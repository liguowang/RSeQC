#!/usr/bin/env python3
"""Element-wise operations on NumPy arrays."""

from __future__ import annotations
import numpy as np

__author__ = "Liguo Wang"

def check_list(v1, v2):
    """Validate that two array-like inputs have the same shape.

    The historical implementation checked only ``.size``. Shape checking is
    safer because arrays such as (2, 3) and (6,) have the same number of
    elements but are not pairwise-compatible in the intended sense.
    """
    a = np.asarray(v1)
    b = np.asarray(v2)

    if a.shape != b.shape:
        raise ValueError(
            f"both arrays must have the same shape: {a.shape} != {b.shape}"
        )

    return a, b


def Add(v1, v2):
    """Return the element-wise sum of two arrays."""
    a, b = check_list(v1, v2)
    return a + b


def Subtract(v1, v2):
    """Return the element-wise difference ``v1 - v2``."""
    a, b = check_list(v1, v2)
    return a - b


def Product(v1, v2):
    """Return the element-wise product of two arrays."""
    a, b = check_list(v1, v2)
    return a * b


def Division(v1, v2):
    """Return element-wise ``(v1 + 1) / (v2 + 1)``.

    Adding 1 to both arrays preserves the behavior documented by the original
    RSeQC implementation.
    """
    a, b = check_list(v1, v2)
    return (a + 1) / (b + 1)


def Average(v1, v2):
    """Return the element-wise arithmetic mean of two arrays."""
    a, b = check_list(v1, v2)
    return (a + b) / 2.0


def geometricMean(v1, v2):
    """Return the element-wise geometric mean of two arrays."""
    a, b = check_list(v1, v2)
    return np.sqrt(a * b)


def Max(v1, v2):
    """Return the element-wise maximum of two arrays."""
    a, b = check_list(v1, v2)
    return np.maximum(a, b)


def Min(v1, v2):
    """Return the element-wise minimum of two arrays."""
    a, b = check_list(v1, v2)
    return np.minimum(a, b)


def euclidean_distance(v1, v2):
    """Return the root-mean-square difference between two arrays.

    Despite its historical name, this function computes
    ``sqrt(mean((v1 - v2)**2))`` rather than the conventional Euclidean
    distance ``sqrt(sum((v1 - v2)**2))``. The historical numerical behavior
    is intentionally preserved.
    """
    a, b = check_list(v1, v2)
    if a.size == 0:
        raise ValueError("arrays must not be empty")
    return np.sqrt(np.mean((a - b) ** 2))
