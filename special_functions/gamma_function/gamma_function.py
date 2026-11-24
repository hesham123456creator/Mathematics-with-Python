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

plt.savefig("special_functions/gamma_function/gamma_function_plot.png", dpi=300)
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
    "special_functions/gamma_function/gamma_function_minimum.png",
    dpi=300
)

plt.show()

print("Improved plot saved successfully.")
