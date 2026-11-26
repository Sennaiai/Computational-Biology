"""Input contracts shared by this study's independent model components."""

import math

import numpy as np


def number(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(f"{name} must be a finite real number")
    try:
        valid = math.isfinite(value) and (value > 0 if positive else value >= 0)
    except (TypeError, ValueError):
        valid = False
    if not valid:
        raise ValueError(
            f"{name} must be finite and {'positive' if positive else 'nonnegative'}"
        )


def integer(value, name, low, high):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer")
    if not low <= value <= high:
        raise ValueError(f"{name} must be between {low} and {high}")


def time_grid(times):
    times = np.asarray(times, dtype=float)
    if (
        times.ndim != 1
        or times.size < 2
        or not np.isfinite(times).all()
        or times[0] != 0
        or np.any(np.diff(times) <= 0)
    ):
        raise ValueError("times must be finite, start at zero and increase strictly")
    return times
