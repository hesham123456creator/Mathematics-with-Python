"""
Lerch Transcendent
==================

Definition:
    Phi(z, s, a) = sum_{n=0}^infinity z^n / (n+a)^s

This module provides a numerical implementation
of the Lerch transcendent using a truncated series.
"""

import numpy as np


def lerch_phi(z, s, a, n_terms=10000):
    """
    Compute the Lerch transcendent using a truncated series.

    Parameters
    ----------
    z : float
        Series variable.
    s : float
        Order parameter.
    a : float
        Shift parameter (must be positive).
    n_terms : int
        Number of terms in the series.

    Returns
    -------
    float
        Numerical approximation of Phi(z, s, a).
    """

    if a <= 0:
        raise ValueError("Parameter a must be positive.")

    if n_terms <= 0:
        raise ValueError("n_terms must be positive.")

    if abs(z) >= 1:
        raise ValueError("This implementation requires |z| < 1.")

    n = np.arange(n_terms)

    terms = z**n / (n + a)**s

    return np.sum(terms)


if __name__ == "__main__":

    z = 0.5
    s = 2
    a = 1

    result = lerch_phi(z, s, a)

    print("Lerch Transcendent")
    print("------------------")
    print(f"z = {z}, s = {s}, a = {a}")
    print(f"Phi(z, s, a) = {result:.12f}")
