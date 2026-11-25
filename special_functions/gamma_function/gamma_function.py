"""
=========================================================
Gamma Function: Mathematical Properties and Computation
=========================================================

Definition:
    Gamma(x) = integral from 0 to infinity of
               t**(x-1) * exp(-t) dt

Fundamental properties:

1. Recurrence relation:
       Gamma(x+1) = x * Gamma(x)

2. Factorial relation:
       Gamma(n) = (n-1)!

3. Special value:
       Gamma(1/2) = sqrt(pi)

Libraries:
    NumPy
    SciPy
    SymPy
    Matplotlib
"""
from scipy.optimize import minimize_scalar
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
import sympy as sp

# ---------------------------------------------------------
# 1. Numerical Evaluation of the Gamma Function
# ---------------------------------------------------------

print("=" * 60)
print("NUMERICAL EVALUATION OF THE GAMMA FUNCTION")
print("=" * 60)

x_values = [0.5, 1, 2, 3, 4, 5, 6]

print(f"{'x':>10} {'Gamma(x)':>20}")

for x in x_values:
    result = gamma(x)
    print(f"{x:10.2f} {result:20.10f}")

# ---------------------------------------------------------
# 2. Verification of the Recurrence Relation
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("VERIFICATION OF THE RECURRENCE RELATION")
print("=" * 60)

test_values = [0.5, 1, 2, 3, 4, 5]

print(f"{'x':>10} {'Gamma(x+1)':>20} {'x*Gamma(x)':>20}")

for x in test_values:
    left_side = gamma(x + 1)
    right_side = x * gamma(x)

    print(f"{x:10.2f} {left_side:20.10f} {right_side:20.10f}")

# ---------------------------------------------------------
# 3. Verification of the Factorial Relation
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("VERIFICATION OF THE FACTORIAL RELATION")
print("=" * 60)

integer_values = [1, 2, 3, 4, 5, 6, 7]

print(f"{'n':>10} {'Gamma(n)':>20} {'(n-1)!':>20}")

for n in integer_values:
    gamma_value = gamma(n)
    factorial_value = math.factorial(n - 1)

    print(f"{n:10d} {gamma_value:20.10f} {factorial_value:20d}")

# ---------------------------------------------------------
# 4. Verification of the Special Value Gamma(1/2)
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("VERIFICATION OF THE SPECIAL VALUE")
print("=" * 60)

gamma_half = gamma(0.5)
sqrt_pi = math.sqrt(math.pi)

print(f"Gamma(1/2) = {gamma_half:.12f}")
print(f"sqrt(pi)   = {sqrt_pi:.12f}")
print(f"Difference = {abs(gamma_half - sqrt_pi):.12e}")

# ---------------------------------------------------------
# 5. Plot of the Gamma Function
# ---------------------------------------------------------

print("\nGenerating the Gamma function plot...")

x_plot = np.linspace(0.1, 5, 500)
y_plot = gamma(x_plot)

plt.figure(figsize=(10, 6))
plt.plot(x_plot, y_plot, label=r"$\Gamma(x)$", linewidth=2)

plt.xlabel("x", fontsize=12)
plt.ylabel(r"$\Gamma(x)$", fontsize=12)
plt.title("The Gamma Function", fontsize=14)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("gamma_function_plot.png", dpi=300)
plt.show()

print("Plot saved successfully.")

# ---------------------------------------------------------
# 6. Numerical Determination of the Minimum
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("NUMERICAL MINIMUM OF THE GAMMA FUNCTION")
print("=" * 60)

result = minimize_scalar(gamma, bounds=(1, 2), method="bounded")

x_min = result.x
gamma_min = result.fun

print(f"Minimum location: x = {x_min:.10f}")
print(f"Minimum value:    Gamma(x) = {gamma_min:.10f}")

# ---------------------------------------------------------
# 7. Improved Gamma Function Plot with Minimum Point
# ---------------------------------------------------------

plt.figure(figsize=(10, 6))

plt.plot(x_plot, y_plot, label=r"$\Gamma(x)$", linewidth=2)

plt.scatter(x_min, gamma_min, color="red", s=70,
            zorder=5, label="Minimum point")

plt.axvline(x_min, color="red", linestyle="--", alpha=0.6)
plt.axhline(gamma_min, color="red", linestyle="--", alpha=0.6)

plt.annotate(
    f"Minimum\n({x_min:.4f}, {gamma_min:.4f})",
    xy=(x_min, gamma_min),
    xytext=(2.0, 1.5),
    arrowprops=dict(arrowstyle="->", color="red"),
    fontsize=10
)

plt.xlabel("x", fontsize=12)
plt.ylabel(r"$\Gamma(x)$", fontsize=12)
plt.title("Gamma Function and Its Numerical Minimum", fontsize=14)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(
    "gamma_function_minimum.png",
    dpi=300
)

plt.show()

print("Improved plot saved successfully.")


# ---------------------------------------------------------
# 8. Euler's Reflection Formula
#
# Gamma(x) * Gamma(1-x) = pi / sin(pi*x)
#
# Valid for 0 < x < 1
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("VERIFICATION OF EULER'S REFLECTION FORMULA")
print("=" * 60)


def reflection_formula(x):
    """Verify Euler's reflection formula for 0 < x < 1."""
    if not 0 < x < 1:
        raise ValueError("x must satisfy 0 < x < 1")

    left_side = gamma(x) * gamma(1 - x)
    right_side = math.pi / math.sin(math.pi * x)

    return left_side, right_side


reflection_values = [0.1, 0.2, 0.3, 0.4, 0.5,
                     0.6, 0.7, 0.8, 0.9]

print(f"{'x':>10} {'Left side':>18} "
      f"{'Right side':>18} {'Absolute error':>18}")

for x in reflection_values:
    left, right = reflection_formula(x)
    error = abs(left - right)

    print(f"{x:10.2f} {left:18.10f} "
          f"{right:18.10f} {error:18.4e}")


# ---------------------------------------------------------
# 9. Visualization of Euler's Reflection Formula
# ---------------------------------------------------------

x_reflection = np.linspace(0.05, 0.95, 300)

left_values = gamma(x_reflection) * gamma(1 - x_reflection)

right_values = np.pi / np.sin(np.pi * x_reflection)

plt.figure(figsize=(10, 6))

plt.plot(x_reflection, left_values,
         label=r"$\Gamma(x)\Gamma(1-x)$", linewidth=2)

plt.plot(x_reflection, right_values, "--",
         label=r"$\pi/\sin(\pi x)$", linewidth=2)

plt.xlabel("x")
plt.ylabel("Function value")
plt.title("Euler's Reflection Formula")

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("gamma_reflection_formula.png", dpi=300)
plt.show()

print("Reflection formula plot saved successfully.")


# ---------------------------------------------------------
# 10. Legendre's Duplication Formula
#
# Gamma(x) * Gamma(x + 1/2)
# = 2^(1-2x) * sqrt(pi) * Gamma(2x)
#
# Valid for x > 0
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("VERIFICATION OF LEGENDRE'S DUPLICATION FORMULA")
print("=" * 60)


def duplication_formula(x):
    """Verify Legendre's duplication formula for x > 0."""
    if x <= 0:
        raise ValueError("x must be positive")

    left_side = gamma(x) * gamma(x + 0.5)

    right_side = (
        2**(1 - 2*x)
        * math.sqrt(math.pi)
        * gamma(2*x)
    )

    return left_side, right_side


duplication_values = [
    0.25, 0.5, 0.75, 1.0, 1.25,
    1.5, 1.75, 2.0, 2.5, 3.0
]

print(f"{'x':>10} {'Left side':>18} "
      f"{'Right side':>18} {'Absolute error':>18}")

for x in duplication_values:
    left, right = duplication_formula(x)
    error = abs(left - right)

    print(f"{x:10.2f} {left:18.10f} "
          f"{right:18.10f} {error:18.4e}")


# ---------------------------------------------------------
# 11. Visualization of Legendre's Duplication Formula
# ---------------------------------------------------------

x_duplication = np.linspace(0.1, 3, 300)

left_values = (
    gamma(x_duplication)
    * gamma(x_duplication + 0.5)
)

right_values = (
    2**(1 - 2*x_duplication)
    * math.sqrt(math.pi)
    * gamma(2*x_duplication)
)

plt.figure(figsize=(10, 6))

plt.plot(
    x_duplication,
    left_values,
    label=r"$\Gamma(x)\Gamma(x+1/2)$",
    linewidth=2
)

plt.plot(
    x_duplication,
    right_values,
    "--",
    label=r"$2^{1-2x}\sqrt{\pi}\Gamma(2x)$",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("Function value")
plt.title("Legendre's Duplication Formula")

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("gamma_duplication_formula.png", dpi=300)
plt.show()

print("Duplication formula plot saved successfully.")


# ---------------------------------------------------------
# 12. Numerical Error Analysis
#
# Comparison of Euler's Reflection Formula and
# Legendre's Duplication Formula
#
# Absolute error:
# E_abs = |L - R|
#
# Relative error:
# E_rel = |L - R| / |R|
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("NUMERICAL ERROR ANALYSIS")
print("=" * 60)

# Reflection formula errors

reflection_errors = []
reflection_relative_errors = []

for x in reflection_values:
    left, right = reflection_formula(x)

    absolute_error = abs(left - right)
    relative_error = absolute_error / abs(right)

    reflection_errors.append(absolute_error)
    reflection_relative_errors.append(relative_error)


# Duplication formula errors

duplication_errors = []
duplication_relative_errors = []

for x in duplication_values:
    left, right = duplication_formula(x)

    absolute_error = abs(left - right)
    relative_error = absolute_error / abs(right)

    duplication_errors.append(absolute_error)
    duplication_relative_errors.append(relative_error)


# Display maximum errors

print("\nREFLECTION FORMULA")
print("-" * 60)

print(
    "Maximum absolute error:",
    f"{max(reflection_errors):.12e}"
)

print(
    "Maximum relative error:",
    f"{max(reflection_relative_errors):.12e}"
)


print("\nDUPLICATION FORMULA")
print("-" * 60)

print(
    "Maximum absolute error:",
    f"{max(duplication_errors):.12e}"
)

print(
    "Maximum relative error:",
    f"{max(duplication_relative_errors):.12e}"
)


# ---------------------------------------------------------
# 13. Comparison of Relative Errors
# ---------------------------------------------------------

plt.figure(figsize=(10, 6))

plt.plot(
    reflection_values,
    reflection_relative_errors,
    marker="o",
    linewidth=2,
    label="Reflection formula"
)

plt.plot(
    duplication_values,
    duplication_relative_errors,
    marker="s",
    linewidth=2,
    label="Duplication formula"
)

plt.xlabel("x", fontsize=12)
plt.ylabel("Relative error", fontsize=12)

plt.title(
    "Numerical Accuracy of Gamma Function Identities",
    fontsize=14
)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(
    "gamma_formula_error_comparison.png",
    dpi=300
)

plt.show()

print("Error comparison plot saved successfully.")
