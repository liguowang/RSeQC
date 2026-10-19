#!/usr/bin/env python3
"""Small statistical utility functions used by RSeQC."""

from __future__ import annotations

import math


__author__ = "Liguo Wang"


def _positive_values(arg):
    """Parse positive numeric values from a comma-separated string."""
    values = []
    for item in arg.split(","):
        item = item.strip()
        if not item:
            continue
        value = float(item)
        if value > 0:
            values.append(value)
    return values


def RSS(arg):
    """Return the root sum of squares of comma-separated integers."""
    values = [int(item.strip()) for item in arg.split(",") if item.strip()]
    return math.sqrt(sum(value * value for value in values))


def H_mean(arg):
    """Return the harmonic mean of non-zero comma-separated values.

    Zero values are ignored, preserving historical RSeQC behavior.
    ``"NA"`` is returned if no non-zero values remain.
    """
    values = [
        float(item.strip())
        for item in arg.split(",")
        if item.strip() and float(item.strip()) != 0
    ]
    if not values:
        return "NA"
    return len(values) / sum(1.0 / value for value in values)


def shannon_entropy(arg):
    """Calculate Shannon entropy (Shannon-Wiener index).

    ``arg`` is an iterable of numeric values. Non-positive values are ignored,
    matching the historical implementation.
    """
    values = [float(value) for value in arg if float(value) > 0]
    if not values:
        return 0

    total = sum(values)
    entropy = sum(
        (value / total) * math.log(value / total)
        for value in values
    )
    return -entropy if entropy != 0 else 0


def shannon_entropy_es(arg):
    """Calculate the Chao-Shen (2003) estimator of Shannon entropy.

    The historical RSeQC algorithm is preserved, including its appended
    pseudo-count of 2.
    """
    values = [float(value) for value in arg if float(value) > 0]
    if not values:
        return "NA"
    if len(values) == 1:
        return 0

    values.append(2)
    singleton = sum(value == 1 for value in values)
    total = sum(values)
    c_bar = 1 - singleton / total

    entropy = 0.0
    for value in values:
        probability = c_bar * value / total
        denominator = 1 - (1 - probability) ** total
        if denominator == 0:
            continue
        entropy += probability * math.log(probability) / denominator

    return -entropy if entropy != 0 else 0


def shannon_entropy_ht(arg):
    """Calculate a Horvitz-Thompson-based Shannon entropy estimator."""
    values = _positive_values(arg)
    if not values:
        return "NA"
    if len(values) == 1:
        return 0

    total = sum(values)
    entropy = 0.0
    for value in values:
        probability = value / total
        denominator = 1 - (1 - probability) ** total
        if denominator == 0:
            continue
        entropy += probability * math.log(probability) / denominator

    return -entropy


def simpson_index(arg):
    """Calculate the Gini-Simpson index from comma-separated values."""
    values = _positive_values(arg)
    if not values:
        return 0

    total = sum(values)
    return 1 - sum((value / total) ** 2 for value in values)


def simpson_index_es(arg):
    """Calculate the unbiased estimator of the Gini-Simpson index."""
    values = _positive_values(arg)
    total = sum(values)

    if total <= 1:
        return 0

    numerator = sum(value * (value - 1) for value in values)
    return 1 - numerator / (total * (total - 1))


def Hill_number(arg, qvalue=1):
    """Calculate Hill diversity for comma-separated abundance values.

    For q=1, return exp(Shannon entropy). For q != 1, use the standard
    Hill-number expression.
    """
    values = _positive_values(arg)
    if not values:
        return 0

    if qvalue == 1:
        return math.exp(shannon_entropy(values))

    total = sum(values)
    frequency_power_sum = sum(
        (value / total) ** qvalue for value in values
    )
    return frequency_power_sum ** (1.0 / (1.0 - qvalue))


def percentile(N, percent, key=lambda x: x):
    """Return a linearly interpolated percentile from an already sorted list."""
    if not N:
        return None
    if not 0 <= percent <= 100:
        raise ValueError("percent must be between 0 and 100")

    k = (len(N) - 1) * percent / 100.0
    lower = math.floor(k)
    upper = math.ceil(k)

    if lower == upper:
        return key(N[int(k)])

    d0 = key(N[lower]) * (upper - k)
    d1 = key(N[upper]) * (k - lower)
    return d0 + d1


def percentile_list(N):
    """Return percentile positions from an already sorted sequence.

    Historical behavior is preserved: empty input returns ``None``; input
    shorter than 100 values is returned unchanged; otherwise percentiles
    1 through 100 are returned as rounded integers.
    """
    if not N:
        return None
    if len(N) < 100:
        return N

    per_list = []
    for percent in range(1, 101):
        k = (len(N) - 1) * percent / 100.0
        lower = math.floor(k)
        upper = math.ceil(k)

        if lower == upper:
            per_list.append(int(N[int(k)]))
        else:
            d0 = N[lower] * (upper - k)
            d1 = N[upper] * (k - lower)
            per_list.append(int(round(d0 + d1)))

    return per_list
