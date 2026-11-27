"""
Mathematical Property Verification
==================================
Lerch Transcendent
"""

import mpmath as mp

mp.mp.dps = 40


def verify_monotonicity_z():
    print("\n1. Monotonicity with respect to z")

    s, a = 2, 1
    z_values = [-0.8, -0.4, 0, 0.4, 0.8]

    values = [mp.lerchphi(z, s, a) for z in z_values]

    for z, value in zip(z_values, values):
        print(f"z = {z:5.2f}, Phi = {mp.nstr(value, 15)}")

    assert all(values[i] < values[i + 1]
               for i in range(len(values) - 1))

    print("Strictly increasing: VERIFIED")


def verify_monotonicity_s():
    print("\n2. Monotonicity with respect to s")

    z, a = mp.mpf("0.5"), 1
    s_values = [0.5, 1, 2, 3, 4, 5]

    values = [mp.lerchphi(z, s, a) for s in s_values]

    for s, value in zip(s_values, values):
        print(f"s = {s:4.1f}, Phi = {mp.nstr(value, 15)}")

    assert all(values[i] > values[i + 1]
               for i in range(len(values) - 1))

    print("Strictly decreasing: VERIFIED")


def verify_derivative_a():
    print("\n3. Derivative identity with respect to a")

    z = mp.mpf("0.5")
    s = mp.mpf("2")
    a = mp.mpf("1")
    h = mp.mpf("1e-20")

    numerical = (
        mp.lerchphi(z, s, a + h)
        - mp.lerchphi(z, s, a - h)
    ) / (2 * h)

    theoretical = -s * mp.lerchphi(z, s + 1, a)

    error = abs(numerical - theoretical)

    print("Numerical derivative =", mp.nstr(numerical, 25))
    print("Theoretical derivative =", mp.nstr(theoretical, 25))
    print("Absolute difference =", mp.nstr(error, 10))

    assert error < mp.mpf("1e-15")

    print("Derivative identity: VERIFIED")


if __name__ == "__main__":

    print("=" * 55)
    print("LERCH TRANSCENDENT: MATHEMATICAL PROPERTY VERIFICATION")
    print("=" * 55)

    verify_monotonicity_z()
    verify_monotonicity_s()
    verify_derivative_a()

    print("\nAll tests completed successfully.")
