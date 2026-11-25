
"""
=========================================================
Digamma Function: Mathematical Properties and Computation
=========================================================

Definition:
    psi(x) = Gamma'(x) / Gamma(x)

The Digamma function is the logarithmic derivative
of the Gamma function.

Fundamental property:
    psi(x+1) = psi(x) + 1/x

Libraries:
    NumPy
    SciPy
    SymPy
    Matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import math

from scipy.special import digamma


# ---------------------------------------------------------
# 1. Numerical Evaluation of the Digamma Function
# ---------------------------------------------------------

print("=" * 60)
print("NUMERICAL EVALUATION OF THE DIGAMMA FUNCTION")
print("=" * 60)

x_values = [0.5, 1, 1.5, 2, 3, 4, 5, 6]

print(f"{'x':>10} {'psi(x)':>20}")

for x in x_values:
    result = digamma(x)
    print(f"{x:10.2f} {result:20.10f}")


# ---------------------------------------------------------
# 2. Verification of the Recurrence Relation
#
# psi(x+1) = psi(x) + 1/x
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("VERIFICATION OF THE RECURRENCE RELATION")
print("=" * 60)

test_values = [0.5, 1, 1.5, 2, 3, 4, 5]

print(f"{'x':>10} {'psi(x+1)':>20} "
      f"{'psi(x)+1/x':>20}")

for x in test_values:
    left_side = digamma(x + 1)
    right_side = digamma(x) + 1 / x

    print(f"{x:10.2f} {left_side:20.10f} "
          f"{right_side:20.10f}")


# ---------------------------------------------------------
# 3. Verification of Special Values
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("VERIFICATION OF SPECIAL VALUES")
print("=" * 60)

# Euler-Mascheroni constant
euler_gamma = -digamma(1.0)

# First special value: psi(1) = -gamma

psi_one = digamma(1.0)

print("\n1. Verification of psi(1) = -gamma")

print(f"psi(1)       = {psi_one:.12f}")
print(f"-gamma       = {-euler_gamma:.12f}")
print(f"Difference   = {abs(psi_one + euler_gamma):.4e}")


# Second special value:
# psi(1/2) = -gamma - 2*ln(2)

psi_half = digamma(0.5)

theoretical_value = -euler_gamma - 2 * math.log(2)

print("\n2. Verification of psi(1/2) = -gamma - 2*ln(2)")

print(f"psi(1/2)     = {psi_half:.12f}")
print(f"Theoretical  = {theoretical_value:.12f}")
print(f"Difference   = {abs(psi_half - theoretical_value):.4e}")


# ---------------------------------------------------------
# 4. Graphical Representation of the Digamma Function
# ---------------------------------------------------------

print("\nGenerating the Digamma function plot...")

x_plot = np.linspace(0.1, 6, 500)

y_plot = digamma(x_plot)

plt.figure(figsize=(10, 6))

plt.plot(
    x_plot,
    y_plot,
    label=r"$\psi(x)$",
    linewidth=2
)

plt.axhline(0, color="black", linewidth=0.8)

plt.xlabel("x", fontsize=12)
plt.ylabel(r"$\psi(x)$", fontsize=12)

plt.title("The Digamma Function", fontsize=14)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("digamma_function_plot.png", dpi=300)
plt.show()

print("Digamma plot saved successfully.")


# ---------------------------------------------------------
# 5. Trigamma Function
#
# psi_1(x) = psi'(x)
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("NUMERICAL EVALUATION OF THE TRIGAMMA FUNCTION")
print("=" * 60)

from scipy.special import polygamma

x_values = [0.5, 1, 1.5, 2, 3, 4, 5, 6]

print(f"{'x':>10} {'psi_1(x)':>20}")

for x in x_values:
    result = polygamma(1, x)
    print(f"{x:10.2f} {result:20.10f}")


# ---------------------------------------------------------
# 6. Graphical Representation of the Trigamma Function
# ---------------------------------------------------------

print("\nGenerating the Trigamma function plot...")

trigamma_values = polygamma(1, x_plot)

plt.figure(figsize=(10, 6))

plt.plot(
    x_plot,
    trigamma_values,
    label=r"$\psi_1(x)$",
    linewidth=2
)

plt.xlabel("x", fontsize=12)
plt.ylabel(r"$\psi_1(x)$", fontsize=12)

plt.title("The Trigamma Function", fontsize=14)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("trigamma_function_plot.png", dpi=300)
plt.show()

print("Trigamma plot saved successfully.")


# ---------------------------------------------------------
# 7. Verification of the Trigamma Recurrence Relation
#
# psi_1(x+1) = psi_1(x) - 1/x^2
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("VERIFICATION OF THE TRIGAMMA RECURRENCE RELATION")
print("=" * 60)

test_values = [0.5, 1, 1.5, 2, 3, 4, 5]

print(f"{'x':>10} {'psi_1(x+1)':>20} "
      f"{'psi_1(x)-1/x^2':>22} {'Error':>15}")

for x in test_values:
    left_side = polygamma(1, x + 1)
    right_side = polygamma(1, x) - 1 / x**2

    error = abs(left_side - right_side)

    print(f"{x:10.2f} {left_side:20.10f} "
          f"{right_side:22.10f} {error:15.4e}")

# ---------------------------------------------------------
# 8. Complete Monotonicity of the Trigamma Function
#
# (-1)^n * psi_1^(n)(x) > 0
#
# psi_1^(n)(x) = polygamma(n+1, x)
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("COMPLETE MONOTONICITY OF THE TRIGAMMA FUNCTION")
print("=" * 60)

orders = range(5)
test_values = [0.5, 1.0, 2.0, 3.0, 5.0]

for n in orders:

    print(f"\nDerivative order n = {n}")
    print("-" * 60)

    print(f"{'x':>10} {'(-1)^n psi_1^(n)(x)':>25}")

    for x in test_values:

        derivative = polygamma(n + 1, x)

        signed_derivative = (-1)**n * derivative

        print(f"{x:10.2f} {signed_derivative:25.12f}")

        assert signed_derivative > 0

print("\nAll tested values satisfy complete monotonicity.")

# ---------------------------------------------------------
# 9. Comparison of Series Approximation Methods
#
# Direct partial sum
# Integral-tail correction
# Euler-Maclaurin correction
# ---------------------------------------------------------

print("\n" + "=" * 75)
print("COMPARISON OF TRIGAMMA DERIVATIVE SERIES METHODS")
print("=" * 75)

test_values = [0.5, 1.0, 2.0, 3.0]
N = 100000

for n in range(5):

    p = n + 2

    print(f"\nDerivative order n = {n}")
    print("-" * 75)

    print(
        f"{'x':>6} {'Direct error':>16} "
        f"{'Integral error':>18} {'Euler-Maclaurin error':>24}"
    )

    for x in test_values:

        scipy_value = (-1)**n * polygamma(n + 1, x)

        k = np.arange(N, dtype=float)

        partial_sum = np.sum(1.0 / (x + k)**p)

        factor = math.factorial(n + 1)

        # Method 1: Direct partial sum
        direct_value = factor * partial_sum

        # Method 2: Integral-tail correction
        integral_tail = 1.0 / (
            (p - 1) * (x + N)**(p - 1)
        )

        integral_value = factor * (
            partial_sum + integral_tail
        )

        # Method 3: Euler-Maclaurin correction
        a = x + N

        em_tail = (
            a**(1 - p) / (p - 1)
            + 0.5 * a**(-p)
            + p / 12.0 * a**(-p - 1)
        )

        em_value = factor * (partial_sum + em_tail)

        direct_error = abs(scipy_value - direct_value)
        integral_error = abs(scipy_value - integral_value)
        em_error = abs(scipy_value - em_value)

        print(
            f"{x:6.2f} "
            f"{direct_error:16.4e} "
            f"{integral_error:18.4e} "
            f"{em_error:24.4e}"
        )

print("\nAll three methods compared successfully.")

# ---------------------------------------------------------
# 10. Error Convergence Comparison
# ---------------------------------------------------------

print("\n" + "=" * 65)
print("ERROR CONVERGENCE COMPARISON")
print("=" * 65)

x = 0.5
n = 0
p = n + 2

exact_value = (-1)**n * polygamma(n + 1, x)

N_values = [10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000]

direct_errors = []
integral_errors = []
em_errors = []

for N in N_values:

    k = np.arange(N, dtype=float)

    partial_sum = np.sum(1.0 / (x + k)**p)

    factor = math.factorial(n + 1)

    # Method 1: Direct partial sum
    direct_value = factor * partial_sum

    # Method 2: Integral-tail correction
    a = x + N

    integral_tail = a**(1 - p) / (p - 1)

    integral_value = factor * (partial_sum + integral_tail)

    # Method 3: Euler-Maclaurin correction
    em_tail = (
        a**(1 - p) / (p - 1)
        + 0.5 * a**(-p)
        + p / 12.0 * a**(-p - 1)
    )

    em_value = factor * (partial_sum + em_tail)

    direct_errors.append(abs(exact_value - direct_value))
    integral_errors.append(abs(exact_value - integral_value))
    em_errors.append(abs(exact_value - em_value))

    print(
        f"N = {N:7d} | "
        f"Direct = {direct_errors[-1]:.4e} | "
        f"Integral = {integral_errors[-1]:.4e} | "
        f"Euler-Maclaurin = {em_errors[-1]:.4e}"
    )

# Plot convergence errors
plt.figure(figsize=(10, 6))

plt.loglog(N_values, direct_errors, 'o-', label="Direct partial sum")
plt.loglog(N_values, integral_errors, 's-', label="Integral correction")
plt.loglog(N_values, em_errors, '^-', label="Euler-Maclaurin correction")

plt.xlabel("Number of terms N")
plt.ylabel("Absolute error")
plt.title("Convergence of Trigamma Series Approximations")
plt.grid(True, which="both", linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()

plt.savefig("trigamma_series_convergence.png", dpi=300)
plt.show()

print("\nConvergence plot saved successfully.")
