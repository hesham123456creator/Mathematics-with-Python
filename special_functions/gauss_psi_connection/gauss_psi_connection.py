import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hyp2f1, digamma
from scipy.optimize import approx_fprime


# --------------------------------------------------
# Connection between Gauss hypergeometric and Psi
#
# D(z) = d/da 2F1(a,2;3;z) at a = 1
# --------------------------------------------------


def numerical_derivative(z, h=1e-5):
    """Numerical differentiation with respect to parameter a."""
    return (
        hyp2f1(1 + h, 2, 3, z)
        - hyp2f1(1 - h, 2, 3, z)
    ) / (2 * h)


def harmonic_series(z, N=1000):
    """Harmonic-number series representation."""
    n = np.arange(1, N + 1)
    harmonic = np.cumsum(1 / n)
    terms = 2 * harmonic / (n + 2) * z**n
    return np.sum(terms)


def closed_form(z):
    """Analytical logarithmic representation."""
    return (2 / z**2) * (
        -z
        + 0.5 * np.log(1 - z)**2
        - (1 - z) * np.log(1 - z)
    )


# --------------------------------------------------
# Numerical comparison
# --------------------------------------------------

z_values = [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]

print("\nComparison of three representations")
print("-" * 75)
print(f"{'z':>8} {'Numerical':>18} {'Series':>18} {'Closed form':>18}")

for z in z_values:
    numerical = numerical_derivative(z)
    series = harmonic_series(z)
    analytical = closed_form(z)

    print(
        f"{z:8.2f} {numerical:18.12f} "
        f"{series:18.12f} {analytical:18.12f}"
    )


# --------------------------------------------------
# Absolute error analysis
# --------------------------------------------------

errors = []

for z in z_values:
    numerical = numerical_derivative(z)
    analytical = closed_form(z)
    errors.append(abs(numerical - analytical))

print("\nMaximum absolute error:")
print(f"{max(errors):.6e}")


# --------------------------------------------------
# Graphical comparison
# --------------------------------------------------

z_plot = np.linspace(0.01, 0.9, 200)

series_values = [harmonic_series(z) for z in z_plot]
closed_values = [closed_form(z) for z in z_plot]

plt.figure(figsize=(9, 6))

plt.plot(z_plot, series_values, label="Harmonic series")
plt.plot(
    z_plot,
    closed_values,
    "--",
    label="Closed-form expression"
)

plt.xlabel("z")
plt.ylabel("D(z)")
plt.title("Gauss Hypergeometric–Psi Connection")
plt.legend()
plt.grid(True)

plt.savefig("images/gauss_psi_comparison.png", dpi=300)
plt.show()

# --------------------------------------------------
# Truncation error analysis
# --------------------------------------------------

N_values = [10, 50, 100, 500]
z_plot = np.linspace(0.01, 0.95, 200)

plt.figure(figsize=(10, 6))

for N in N_values:
    errors = []

    for z in z_plot:
        exact = closed_form(z)
        approximation = harmonic_series(z, N)
        errors.append(abs(exact - approximation))

    plt.semilogy(
        z_plot,
        errors,
        label=f"N = {N}"
    )

plt.xlabel("z")
plt.ylabel("Absolute truncation error")
plt.title("Convergence of the Harmonic-Number Series")
plt.legend()
plt.grid(True)

plt.savefig("images/gauss_psi_truncation_error.png", dpi=300)
plt.show()

# --------------------------------------------------
# High-precision asymptotic error verification
# --------------------------------------------------

import mpmath as mp

mp.mp.dps = 800

print("\nHigh-precision asymptotic error verification")
print("-" * 90)
print(
    f"{'z':>8} {'N':>8} "
    f"{'Actual error':>25} {'Estimate':>25} {'Ratio':>15}"
)

for z_value in ["0.3", "0.5", "0.7", "0.9"]:
    z_mp = mp.mpf(z_value)

    for N in [50, 100, 200, 500, 1000]:

        # High-precision harmonic number
        H = mp.digamma(N + 2) + mp.euler

        # Exact function value from the closed form
        exact = (
            2 / z_mp**2
            * (
                -z_mp
                + mp.mpf("0.5") * mp.log(1 - z_mp)**2
                - (1 - z_mp) * mp.log(1 - z_mp)
            )
        )

        # High-precision partial sum
        partial = mp.fsum(
            2 * mp.harmonic(n) / (n + 2) * z_mp**n
            for n in range(1, N + 1)
        )

        actual_error = exact - partial

        # Asymptotic estimate
        estimate = (
            2 * H / (N + 3)
            * z_mp**(N + 1) / (1 - z_mp)
        )

        ratio = actual_error / estimate

        print(
            f"{z_value:>8} {N:8d} "
            f"{mp.nstr(actual_error, 12):>25} "
            f"{mp.nstr(estimate, 12):>25} "
            f"{mp.nstr(ratio, 12):>15}"
        )

# --------------------------------------------------
# Graphical verification of asymptotic behavior
# --------------------------------------------------

mp.mp.dps = 200

z_values = ["0.3", "0.5", "0.7", "0.9"]
N_values = range(10, 301, 10)

plt.figure(figsize=(10, 6))

for z_value in z_values:
    z = mp.mpf(z_value)
    ratios = []

    for N in N_values:

        # Exact value using high precision
        exact = (
            2 / z**2
            * (
                -z
                + mp.mpf("0.5") * mp.log(1 - z)**2
                - (1 - z) * mp.log(1 - z)
            )
        )

        # Partial sum
        partial = mp.fsum(
            2 * mp.harmonic(n) / (n + 2) * z**n
            for n in range(1, N + 1)
        )

        # Actual remainder
        remainder = exact - partial

        # Asymptotic estimate
        H = mp.harmonic(N + 1)

        estimate = (
            2 * H / (N + 3)
            * z**(N + 1) / (1 - z)
        )

        ratios.append(float(remainder / estimate))

    plt.plot(
        list(N_values),
        ratios,
        marker="o",
        markersize=3,
        label=f"z = {z_value}"
    )

plt.axhline(
    y=1,
    linestyle="--",
    label="Theoretical limit = 1"
)

plt.xlabel("Number of terms N")
plt.ylabel("Actual remainder / Asymptotic estimate")
plt.title("Asymptotic Convergence of the Harmonic-Number Series")
plt.legend()
plt.grid(True)

plt.savefig(
   "images/gauss_psi_asymptotic_ratio.png",    
   dpi=300
)

plt.show()

# --------------------------------------------------
# Comparison of leading and second-order estimates
# --------------------------------------------------

mp.mp.dps = 100

z_values = ["0.3", "0.5", "0.7", "0.9"]
N_values = [50, 100, 200, 500]

print("\nComparison of asymptotic approximations")
print("-" * 100)
print(
    f"{'z':>6} {'N':>6} "
    f"{'Leading error':>18} {'Second-order error':>22}"
)

for z_value in z_values:
    z = mp.mpf(z_value)

    for N in N_values:
        m = N + 1

        a = lambda n: 2 * mp.harmonic(n) / (n + 2)

        am = a(m)

        delta = a(m + 1) / am - 1

        epsilon = (
            a(m + 2) - 2 * a(m + 1) + am
        ) / am

        exact = (
            2 / z**2
            * (
                -z
                + mp.mpf("0.5") * mp.log(1 - z)**2
                - (1 - z) * mp.log(1 - z)
            )
        )

        partial = mp.fsum(
            a(n) * z**n
            for n in range(1, N + 1)
        )

        actual = mp.fsum(
           a(n) * z**n
           for n in range(N + 1, N + 3001)
        )

        leading = am * z**m / (1 - z)

        second_order = am * z**m * (
            1 / (1 - z)
            + delta * z / (1 - z)**2
            + epsilon * z**2 / (1 - z)**3
        )

        leading_error = abs(actual - leading) / abs(actual)

        second_error = (
            abs(actual - second_order) / abs(actual)
        )

        print(
            f"{z_value:>6} {N:6d} "
            f"{float(leading_error):18.6e} "
            f"{float(second_error):22.6e}"
        )


# --------------------------------------------------
# Graphical comparison
# --------------------------------------------------

z = mp.mpf("0.9")
N_plot = list(range(20, 301, 10))

leading_errors = []
second_errors = []

for N in N_plot:
    m = N + 1

    a = lambda n: 2 * mp.harmonic(n) / (n + 2)

    am = a(m)
    delta = a(m + 1) / am - 1
    epsilon = (a(m + 2) - 2 * a(m + 1) + am) / am

    exact = (
        2 / z**2
        * (
            -z
            + mp.mpf("0.5") * mp.log(1 - z)**2
            - (1 - z) * mp.log(1 - z)
        )
    )

    partial = mp.fsum(
        a(n) * z**n for n in range(1, N + 1)
    )

    actual = mp.fsum(
       a(n) * z**n
       for n in range(N + 1, N + 3001)
    )

    leading = am * z**m / (1 - z)

    second_order = am * z**m * (
        1 / (1 - z)
        + delta * z / (1 - z)**2
        + epsilon * z**2 / (1 - z)**3
    )

    leading_errors.append(
        float(abs(actual - leading) / abs(actual))
    )

    second_errors.append(
        float(abs(actual - second_order) / abs(actual))
    )

plt.figure(figsize=(10, 6))

plt.semilogy(
    N_plot,
    leading_errors,
    "o-",
    label="Leading-order approximation"
)

plt.semilogy(
    N_plot,
    second_errors,
    "s-",
    label="Second-order approximation"
)

plt.xlabel("Number of terms N")
plt.ylabel("Relative error")
plt.title("Accuracy of Asymptotic Approximations (z = 0.9)")
plt.legend()
plt.grid(True)

plt.savefig(
    "images/gauss_psi_asymptotic_accuracy.png",
    dpi=300
)

plt.show()
