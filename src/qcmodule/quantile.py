#!/usr/bin/env python3
"""Compute sample quantiles using the nine Hyndman-Fan definitions.

Based on ``quantile.py`` by Ernesto P. Adorio, PhD
(original version 0.0.1, August 7, 2009).
"""

from __future__ import annotations
import math


_QUANTILE_PARAMETERS = (
    (0.0, 0.0, 1.0, 0.0),               # type 1
    (0.5, 0.0, 1.0, 0.0),               # type 2
    (0.5, 0.0, 0.0, 0.0),               # type 3
    (0.0, 0.0, 0.0, 1.0),               # type 4
    (0.5, 0.0, 0.0, 1.0),               # type 5
    (0.0, 1.0, 0.0, 1.0),               # type 6
    (1.0, -1.0, 0.0, 1.0),              # type 7
    (1.0 / 3.0, 1.0 / 3.0, 0.0, 1.0),  # type 8
    (3.0 / 8.0, 0.25, 0.0, 1.0),        # type 9
)


def quantile(x, q, qtype=7, issorted=False):
    """Return a sample quantile.

    Parameters
    ----------
    x : sequence of numbers
        Input observations.
    q : float
        Quantile probability in the closed interval [0, 1].
    qtype : int, default=7
        Hyndman-Fan/R quantile type, from 1 through 9.
    issorted : bool, default=False
        If True, assume ``x`` is already sorted.

    Returns
    -------
    number
        Requested sample quantile.

    Raises
    ------
    TypeError
        If ``qtype`` is not an integer.
    ValueError
        If ``x`` is empty, ``q`` is outside [0, 1], or ``qtype`` is
        outside 1..9.
    """
    if not isinstance(qtype, int) or isinstance(qtype, bool):
        raise TypeError("qtype must be an integer from 1 through 9")
    if not 1 <= qtype <= 9:
        raise ValueError("qtype must be between 1 and 9")
    if not 0 <= q <= 1:
        raise ValueError("q must be between 0 and 1")

    y = x if issorted else sorted(x)
    n = len(y)
    if n == 0:
        raise ValueError("x must contain at least one value")
    if n == 1:
        return y[0]

    a, b, c, d = _QUANTILE_PARAMETERS[qtype - 1]
    fractional, integer = math.modf(a + (n + b) * q - 1)

    if integer < 0:
        return y[0]
    if integer >= n - 1:
        # If the requested position is at/past the last observation, there
        # is no y[j+1] to interpolate with. The sample maximum is the proper
        # boundary value.
        return y[n - 1]

    j = int(math.floor(integer))
    if fractional == 0:
        return y[j]

    return y[j] + (y[j + 1] - y[j]) * (c + d * fractional)


def _test():
    """Run the example supplied with the original module."""
    x = [11.4, 17.3, 21.3, 25.9, 40.1, 50.5, 60.0, 70.0, 75]
    for qtype in range(1, 10):
        print(qtype, quantile(x, 0.35, qtype))


if __name__ == "__main__":
    _test()
