"""Shared utilities for aviary web components."""


def quality_tier(completeness, contamination):
    """Return 'high', 'medium', or 'low' quality tier for a bin."""
    try:
        c, x = float(completeness), float(contamination)
        if c >= 90 and x <= 5:
            return "high"
        if c >= 50 and x <= 10:
            return "medium"
    except (ValueError, TypeError):
        pass
    return "low"
