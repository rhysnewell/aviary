import pytest
import functools

def pytest_addoption(parser):
    parser.addoption("--run-expensive", action="store_true", help="run expensive tests")
    parser.addoption("--run-qsub", action="store_true", help="run tests that require qsub")

def pytest_configure(config):
    config.addinivalue_line("markers", "expensive: mark test as requiring --run-expensive")
    config.addinivalue_line("markers", "qsub: mark test as requiring qsub")

def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-expensive") and config.getoption("--run-qsub"):
        return

    skip_expensive = pytest.mark.skip(reason="need --run-expensive option")
    skip_qsub = pytest.mark.skip(reason="need --run-qsub option")

    for item in items:
        if "expensive" in item.keywords and not config.getoption("--run-expensive"):
            item.add_marker(skip_expensive)
        if "qsub" in item.keywords and not config.getoption("--run-qsub"):
            item.add_marker(skip_qsub)

def skip_unless_expensive(flag_option="--run-expensive"):
    def decorator(test_func):
        @functools.wraps(test_func)
        def wrapper(*args, **kwargs):
            # Find the 'request' fixture in the arguments
            request = kwargs.get('request', None)
            if request is None:
                # If 'request' is not in kwargs, try to find it in positional args
                for arg in args:
                    if hasattr(arg, "config"):  # crude check for pytest's request
                        request = arg
                        break
            if request is None:
                raise RuntimeError("Test must use the 'request' fixture for skip_unless_expensive to work.")

            if not request.config.getoption(flag_option):
                pytest.skip(f"Requires {flag_option} option")
            return test_func(*args, **kwargs)
        return wrapper
    return decorator

def skip_unless_qsub(flag_option="--run-qsub"):
    def decorator(test_func):
        @functools.wraps(test_func)
        def wrapper(*args, **kwargs):
            # Find the 'request' fixture in the arguments
            request = kwargs.get('request', None)
            if request is None:
                # If 'request' is not in kwargs, try to find it in positional args
                for arg in args:
                    if hasattr(arg, "config"):  # crude check for pytest's request
                        request = arg
                        break
            if request is None:
                raise RuntimeError("Test must use the 'request' fixture for skip_unless_qsub to work.")

            if not request.config.getoption(flag_option):
                pytest.skip(f"Requires {flag_option} option")
            return test_func(*args, **kwargs)
        return wrapper
    return decorator