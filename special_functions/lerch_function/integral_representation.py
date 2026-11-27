"""
Integral Representation of the Lerch Transcendent
==================================================
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma


def lerch_phi_integral(z, s, a):

    if a <= 0:
        raise ValueError("Parameter a must be positive.")

    if s <= 0:
        raise ValueError("Parameter s must be positive.")

    if abs(z) >= 1:
        raise ValueError("This implementation requires |z| < 1.")

    def integrand(t):
        numerator = t**(s - 1) * np.exp(-a * t)
        denominator = 1 - z * np.exp(-t)

        return numerator / denominator

    integral, error = quad(
        integrand,
        0,
        np.inf,
        epsabs=1e-12,
        epsrel=1e-12,
        limit=200
    )

    return integral / gamma(s), error / gamma(s)


if __name__ == "__main__":

    z = 0.5
    s = 2
    a = 1

    result, error = lerch_phi_integral(z, s, a)

    print("Lerch Transcendent: Integral Representation")
    print("--------------------------------------------")
    print(f"z = {z}, s = {s}, a = {a}")
    print(f"Integral result = {result:.15f}")
    print(f"Estimated integration error = {error:.6e}")
