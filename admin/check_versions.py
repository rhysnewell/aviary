#!/usr/bin/env python3

"""
Pixi Dependency Version Checker
Checks if dependencies are up to date across all environments
"""

import json
import subprocess
import sys
import re
import argparse
import asyncio
import time
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from packaging import version
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading


@dataclass
class Package:
    name: str
    current_version: str
    latest_version: Optional[str] = None
    is_up_to_date: Optional[bool] = None


class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    BOLD = '\033[1m'
    NC = '\033[0m'  # No Color


# Thread-safe print lock
print_lock = threading.Lock()

# Cache for package version lookups to avoid repeated searches
version_cache = {}
cache_lock = threading.Lock()


def print_status(color: str, message: str) -> None:
    """Print colored status message (thread-safe)."""
    with print_lock:
        print(f"{color}{message}{Colors.NC}")


def print_thread_safe(message: str) -> None:
    """Thread-safe print function."""
    with print_lock:
        print(message)


def run_command(cmd: List[str]) -> Tuple[bool, str]:
    """Run a command and return success status and output."""
    try:
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            check=False
        )
        return result.returncode == 0, result.stdout.strip()
    except Exception as e:
        return False, str(e)


def get_environments() -> List[str]:
    """Get list of pixi environments."""
    print_status(Colors.BLUE, "🔍 Discovering environments...")
    
    success, output = run_command(["pixi", "workspace", "environment", "list"])
    if not success:
        print_status(Colors.RED, f"❌ Failed to get environments: {output}")
        return []
    
    environments = []
    for line in output.split('\n'):
        if line.strip().startswith('- '):
            env_name = line.strip()[2:].split(':')[0]
            environments.append(env_name)
    
    if environments:
        print_status(Colors.BLUE, f"Found environments: {', '.join(environments)}")
    else:
        print_status(Colors.YELLOW, "⚠️  No environments found")
    
    return environments


def get_environment_dependencies(env_name: str) -> List[Package]:
    """Get dependencies for a specific environment."""
    success, output = run_command(["pixi", "list", "-e", env_name, "--json-pretty", "-x"])
    
    if not success:
        print_status(Colors.RED, f"❌ Failed to get dependencies for {env_name}: {output}")
        return []
    
    try:
        deps_data = json.loads(output)
        packages = []
        
        for dep in deps_data:
            if isinstance(dep, dict) and 'name' in dep and 'version' in dep:
                packages.append(Package(
                    name=dep['name'],
                    current_version=dep['version']
                ))
        
        return packages
    
    except json.JSONDecodeError as e:
        print_status(Colors.RED, f"❌ Failed to parse dependencies JSON for {env_name}: {e}")
        return []


def extract_latest_version(package_name: str, search_output: str) -> Optional[str]:
    """Extract the latest version from pixi search output."""
    lines = search_output.split('\n')
    
    # Look for the first line that starts with package_name-version
    pattern = rf"^{re.escape(package_name)}-([^-]+)-"
    
    for line in lines:
        match = re.match(pattern, line.strip())
        if match:
            return match.group(1)
    
    return None


def get_latest_version(package_name: str) -> Optional[str]:
    """Get the latest version of a package from pixi search with caching."""
    # Check cache first
    with cache_lock:
        if package_name in version_cache:
            return version_cache[package_name]
    
    # Add a small delay to avoid overwhelming conda infrastructure
    # Only delay if this is not the first package (to allow some parallelism)
    if len(version_cache) > 0:
        time.sleep(0.05)
    
    success, output = run_command(["pixi", "search", package_name])
    
    if not success:
        result = None
    else:
        result = extract_latest_version(package_name, output)
    
    # Cache the result (even if None to avoid repeated failed lookups)
    with cache_lock:
        version_cache[package_name] = result
    
    return result


def compare_versions(current: str, latest: str) -> bool:
    """Compare two version strings. Returns True if current >= latest."""
    try:
        return version.parse(current) >= version.parse(latest)
    except version.InvalidVersion:
        # Fallback to string comparison if version parsing fails
        return current == latest


def check_package_version(package: Package) -> Tuple[Package, str]:
    """Check if a package is up to date. Returns (package, status_message)."""
    latest = get_latest_version(package.name)
    
    if latest is None:
        return package, f"  Checking {package.name} ({package.current_version})... ⚠️  Could not check latest version"
    
    package.latest_version = latest
    package.is_up_to_date = compare_versions(package.current_version, latest)
    
    if package.is_up_to_date:
        return package, f"  Checking {package.name} ({package.current_version})... ✅ Up to date"
    else:
        return package, f"  Checking {package.name} ({package.current_version})... ❌ Outdated (latest: {latest})"


def check_environment(env_name: str, max_workers: int = 8) -> Tuple[List[Package], int, int]:
    """Check all packages in an environment in parallel. Returns (packages, total_count, outdated_count)."""
    print_status(Colors.YELLOW, f"📦 Checking environment: {env_name}")
    print_thread_safe("=" * 50)
    
    packages = get_environment_dependencies(env_name)
    
    if not packages:
        print_status(Colors.YELLOW, f"⚠️  No explicit dependencies found in environment: {env_name}")
        print_thread_safe("")
        return [], 0, 0
    
    checked_packages = []
    outdated_count = 0
    
    # Limit concurrent pixi search operations to avoid overwhelming the system
    # Use at most 3 threads for I/O bound pixi search operations
    effective_workers = min(3, max_workers, len(packages))
    
    # Check packages in parallel with limited concurrency
    with ThreadPoolExecutor(max_workers=effective_workers) as executor:
        # Submit all package checks
        future_to_package = {executor.submit(check_package_version, package): package for package in packages}
        
        # Collect results as they complete
        for future in as_completed(future_to_package):
            checked_package, status_message = future.result()
            checked_packages.append(checked_package)
            
            # Print the complete status line
            if "⚠️" in status_message:
                print_status(Colors.YELLOW, status_message)
            elif "✅" in status_message:
                print_status(Colors.GREEN, status_message)
            else:  # ❌
                print_status(Colors.RED, status_message)
            
            if checked_package.is_up_to_date is False:
                outdated_count += 1
    
    print_thread_safe("")
    if outdated_count == 0:
        print_status(Colors.GREEN, f"✅ All {len(packages)} packages in '{env_name}' are up to date!")
    else:
        print_status(Colors.RED, f"❌ {outdated_count} out of {len(packages)} packages in '{env_name}' are outdated")
    
    print_thread_safe("")
    return checked_packages, len(packages), outdated_count


def print_summary(total_packages: int, total_outdated: int) -> None:
    """Print summary of all checks."""
    print("=" * 60)
    print_status(Colors.BLUE, "📊 SUMMARY")
    print("=" * 60)
    
    if total_outdated == 0:
        print_status(Colors.GREEN, "🎉 All dependencies are up to date across all environments!")
    else:
        print_status(Colors.RED, f"⚠️  Found {total_outdated} outdated packages out of {total_packages} total packages")
        print()
        print_status(Colors.YELLOW, "💡 To update packages, you can:")
        print("   • Update individual packages: pixi update <package-name>")
        print("   • Update all packages in an environment: pixi update -e <environment>")
        print("   • Update all packages in all environments: pixi update")
    
    print()
    print_status(Colors.BLUE, "✨ Check complete!")


def main() -> None:
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Check dependency versions across pixi environments"
    )
    parser.add_argument(
        "-e", "--environment",
        action="append",
        dest="environments",
        help="Specific environment(s) to check (can be used multiple times). If not specified, all environments will be checked."
    )
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=8,
        help="Number of parallel threads to use (default: 8)"
    )
    # Accept additional positional arguments as environments
    parser.add_argument(
        "additional_environments",
        nargs="*",
        help="Additional environments to check (space-separated)"
    )
    
    args = parser.parse_args()
    
    print_status(Colors.BLUE, "🔍 Checking dependency versions...")
    print_thread_safe("")
    
    # Check if packaging module is available
    try:
        import packaging.version
    except ImportError:
        print_status(Colors.YELLOW, "⚠️  'packaging' module not found. Version comparison may be less accurate.")
        print_status(Colors.YELLOW, "   Install with: pip install packaging")
        print_thread_safe("")
    
    # Combine environments from -e flags and positional arguments
    environments = []
    if args.environments:
        environments.extend(args.environments)
    if args.additional_environments:
        environments.extend(args.additional_environments)
    
    if environments:
        # Use specified environments
        print_status(Colors.BLUE, f"Checking specified environments: {', '.join(environments)}")
    else:
        # Get all environments
        environments = get_environments()
        
        if not environments:
            print_status(Colors.RED, "❌ No environments found or pixi workspace not available")
            sys.exit(1)
        
        print_status(Colors.BLUE, f"Checking all environments: {', '.join(environments)}")
    
    print_thread_safe("")
    
    total_packages = 0
    total_outdated = 0
    all_packages = {}
    
    # Check environments sequentially (to maintain proper separation)
    for env_name in environments:
        try:
            packages, pkg_count, outdated_count = check_environment(env_name, args.jobs)
            all_packages[env_name] = packages
            total_packages += pkg_count
            total_outdated += outdated_count
        except Exception as exc:
            print_status(Colors.RED, f"❌ Environment {env_name} generated an exception: {exc}")
    
    # Print summary
    print_summary(total_packages, total_outdated)
    
    # Exit with error code if any packages are outdated
    sys.exit(1 if total_outdated > 0 else 0)


if __name__ == "__main__":
    main()