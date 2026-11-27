import mpmath as mp
import matplotlib.pyplot as plt
import os

# --------------------------------------------------
# Rigorous bounds for the Gauss-Psi remainder
#
# R_N(z) = sum_{n=N+1}^\infty a_n z^n
# a_n = 2 H_n / (n+2)
#
# L_p(N,z) < R_N(z) < U_p(N,z)
# --------------------------------------------------

mp.mp.dps = 100

os.makedirs("images", exist_ok=True)


def coefficient(n):
    """Compute a_n = 2 H_n / (n+2)."""
    return 2 * mp.harmonic(n) / (n + 2)


def exact_function(z):
    """Closed-form representation of D(z)."""
    return (2 / z**2) * (
        -z
        + mp.log(1 - z)**2 / 2
        - (1 - z) * mp.log(1 - z)
    )


def partial_sum(N, z):
    """Compute the first N terms of D(z)."""
    return mp.fsum(
        coefficient(n) * z**n
        for n in range(1, N + 1)
    )


def actual_remainder(N, z):
    """High-precision evaluation of R_N(z)."""
    return exact_function(z) - partial_sum(N, z)


def lower_bound(N, z, p):
    """Lower bound L_p(N,z)."""
    return mp.fsum(
        coefficient(N + 1 + k) * z**(N + 1 + k)
        for k in range(p)
    )


def upper_bound(N, z, p):
    """Upper bound U_p(N,z)."""
    L = lower_bound(N, z, p)

    return (
        L
        + coefficient(N + p + 1)
        * z**(N + p + 1) / (1 - z)
    )


# --------------------------------------------------
# Numerical verification
# --------------------------------------------------

print("\nGAUSS-PSI REMAINDER BOUNDS")
print("=" * 75)

N_values = [10, 50, 100]
z_values = [0.3, 0.5, 0.7, 0.9]
p_values = [1, 2, 5, 10, 20]

for z in z_values:
    print(f"\n z = {z}")

    for N in N_values:
        R = actual_remainder(N, mp.mpf(str(z)))

        print(f"\n N = {N}")
        print(f"Actual remainder = {mp.nstr(R, 12)}")

        for p in p_values:
            L = lower_bound(N, mp.mpf(str(z)), p)
            U = upper_bound(N, mp.mpf(str(z)), p)

            verified = L < R < U

            width = U - L

            print(
                f"p={p:2d} | "
                f"L={mp.nstr(L, 8):>14} | "
                f"U={mp.nstr(U, 8):>14} | "
                f"Width={mp.nstr(width, 5):>12} | "
                f"Verified={verified}"
            )


# --------------------------------------------------
# Plot: Enclosure width versus p
# --------------------------------------------------

N = 50
z = mp.mpf("0.7")

p_plot = list(range(1, 31))
widths = [
    float(upper_bound(N, z, p) - lower_bound(N, z, p))
    for p in p_plot
]

plt.figure(figsize=(9, 6))

plt.semilogy(p_plot, widths, marker="o")

plt.xlabel("Number of retained terms p")
plt.ylabel("Width of the enclosure")
plt.title("Convergence of Gauss-Psi Remainder Bounds")

plt.grid(True, which="both", linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig(
    "images/gauss_psi_remainder_bounds.png",
    dpi=300
)

plt.show()

print("\nPlot saved successfully!")

# --------------------------------------------------
# Sharpness study: relative errors of the bounds
# --------------------------------------------------

N = 50
z = mp.mpf("0.7")

p_values = list(range(1, 31))

R = actual_remainder(N, z)

lower_errors = []
upper_errors = []

print("\nSHARPNESS STUDY")
print("=" * 65)
print(f"N = {N}, z = {z}")
print(f"Actual remainder = {mp.nstr(R, 12)}")

for p in p_values:
    L = lower_bound(N, z, p)
    U = upper_bound(N, z, p)

    error_L = (R - L) / R
    error_U = (U - R) / R

    lower_errors.append(float(error_L))
    upper_errors.append(float(error_U))

    if p in [1, 2, 5, 10, 20, 30]:
        print(
            f"p = {p:2d} | "
            f"Lower relative error = {mp.nstr(error_L, 8)} | "
            f"Upper relative error = {mp.nstr(error_U, 8)}"
        )

plt.figure(figsize=(9, 6))

plt.semilogy(
    p_values, lower_errors,
    marker="o", label="Lower bound error"
)

plt.semilogy(
    p_values, upper_errors,
    marker="s", label="Upper bound error"
)

plt.xlabel("Number of retained terms p")
plt.ylabel("Relative error")
plt.title("Sharpness of Gauss-Psi Remainder Bounds")

plt.legend()
plt.grid(True, which="both", linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig(
    "images/gauss_psi_bound_sharpness.png",
    dpi=300
)

plt.show()

print("\nSharpness study completed successfully!")

# --------------------------------------------------
# Sharper harmonic-number upper bound
# --------------------------------------------------

def sharper_upper_bound(N, z, p):
    """
    Sharper upper bound for the Gauss-Psi remainder.

    Uses an estimate involving harmonic numbers.
    """

    m = N + p + 1

    Hm = mp.harmonic(m)

    Cm = Hm - mp.mpf(m + 2) / (m + 1)

    Sm = (
        -mp.log(1 - z)
        - sum(z**k / k for k in range(1, m + 2))
    ) / z**(m + 2)

    L = lower_bound(N, z, p)

    tail_bound = z**m * (
        2 / ((m + 1) * (1 - z))
        + 2 * Cm * Sm
    )

    return L + tail_bound

# --------------------------------------------------
# Comparison: original and sharper upper bounds
# --------------------------------------------------

print("\nCOMPARISON OF UPPER BOUNDS")
print("=" * 75)

test_cases = [
    (10, mp.mpf("0.5"), 5),
    (50, mp.mpf("0.7"), 10),
    (100, mp.mpf("0.9"), 20),
]

for N, z, p in test_cases:

    R = actual_remainder(N, z)

    U_original = upper_bound(N, z, p)
    U_sharper = sharper_upper_bound(N, z, p)

    print(f"\nN = {N}, z = {z}, p = {p}")
    print(f"Actual remainder:  {mp.nstr(R, 12)}")
    print(f"Original upper:    {mp.nstr(U_original, 12)}")
    print(f"Sharper upper:     {mp.nstr(U_sharper, 12)}")

    print(
        "Sharper bound is valid:",
        U_sharper > R
    )

    print(
        "Sharper than original:",
        U_sharper < U_original
    )
# --------------------------------------------------
# Systematic verification of the sharper upper bound
# --------------------------------------------------

print("\nSYSTEMATIC VERIFICATION OF THE NEW THEOREM")
print("=" * 75)

test_N = [1, 5, 10, 50, 100]
test_p = [1, 2, 5, 10, 20]
test_z = [
    mp.mpf("0.2"),
    mp.mpf("0.5"),
    mp.mpf("0.7"),
    mp.mpf("0.9"),
    mp.mpf("0.99")
]

total_tests = 0
passed_tests = 0

for N in test_N:
    for p in test_p:
        for z in test_z:

            R = actual_remainder(N, z)
            L = lower_bound(N, z, p)
            U_new = sharper_upper_bound(N, z, p)
            U_old = upper_bound(N, z, p)

            valid = (
                L < R < U_new < U_old
            )

            total_tests += 1

            if valid:
                passed_tests += 1

print(f"Total tests: {total_tests}")
print(f"Passed tests: {passed_tests}")
print(f"Failed tests: {total_tests - passed_tests}")

if passed_tests == total_tests:
    print("\nALL TESTS PASSED!")
else:
    print("\nSOME TESTS FAILED!")
