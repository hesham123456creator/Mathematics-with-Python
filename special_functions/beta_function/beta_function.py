# =========================================================
# Beta Function
#
# B(x,y) = Gamma(x) Gamma(y) / Gamma(x+y)
#
# Numerical study using Python
# =========================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import beta, gamma
import math


# ---------------------------------------------------------
# 1. Numerical Evaluation of the Beta Function
# ---------------------------------------------------------

print("\n" + "=" * 60)
print("NUMERICAL EVALUATION OF THE BETA FUNCTION")
print("=" * 60)

test_values = [
    (0.5, 0.5),
    (1.0, 1.0),
    (1.0, 2.0),
    (2.0, 3.0),
    (3.0, 4.0),
    (0.5, 2.0),
    (1.5, 2.5)
]

print(f"{'x':>10} {'y':>10} {'B(x,y)':>20}")

for x, y in test_values:

    result = beta(x, y)

    print(f"{x:10.2f} {y:10.2f} {result:20.12f}")
# ---------------------------------------------------------
# 2. Integral Representation of the Beta Function
#
# B(x,y) = integral from 0 to 1 of
#          t^(x-1) * (1-t)^(y-1) dt
# ---------------------------------------------------------

from scipy.integrate import quad

print("\n" + "=" * 65)
print("INTEGRAL REPRESENTATION OF THE BETA FUNCTION")
print("=" * 65)

print(
    f"{'x':>8} {'y':>8} "
    f"{'SciPy beta':>18} {'Numerical integral':>22} {'Error':>14}"
)

for x, y in test_values:

    integrand = lambda t: t**(x - 1) * (1 - t)**(y - 1)

    integral_value, estimated_error = quad(
        integrand, 0, 1
    )

    scipy_value = beta(x, y)

    difference = abs(scipy_value - integral_value)

    print(
        f"{x:8.2f} {y:8.2f} "
        f"{scipy_value:18.12f} "
        f"{integral_value:22.12f} "
        f"{difference:14.4e}"
    )

print("\nIntegral representation verified numerically.")

# ---------------------------------------------------------
# 3. Verification of the Symmetry Property
#
# B(x,y) = B(y,x)
# ---------------------------------------------------------

print("\n" + "=" * 65)
print("VERIFICATION OF THE SYMMETRY PROPERTY")
print("=" * 65)

print(
    f"{'x':>8} {'y':>8} "
    f"{'B(x,y)':>18} {'B(y,x)':>18} {'Error':>14}"
)

for x, y in test_values:

    value_xy = beta(x, y)
    value_yx = beta(y, x)

    error = abs(value_xy - value_yx)

    print(
        f"{x:8.2f} {y:8.2f} "
        f"{value_xy:18.12f} "
        f"{value_yx:18.12f} "
        f"{error:14.4e}"
    )

print("\nSymmetry property verified numerically.")

# ---------------------------------------------------------
# 4. Verification of the Recurrence Relations
#
# B(x+1,y) = x/(x+y) * B(x,y)
# B(x,y+1) = y/(x+y) * B(x,y)
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("VERIFICATION OF THE BETA FUNCTION RECURRENCE RELATIONS")
print("=" * 70)

print("\n1. Verification of B(x+1,y) = x/(x+y) * B(x,y)")
print("-" * 70)
print(f"{'x':>8} {'y':>8} {'Left side':>18} {'Right side':>18} {'Error':>14}")

for x, y in test_values:
    lhs = beta(x + 1, y)
    rhs = (x / (x + y)) * beta(x, y)
    error = abs(lhs - rhs)
    print(f"{x:8.2f} {y:8.2f} {lhs:18.12f} {rhs:18.12f} {error:14.4e}")

print("\n2. Verification of B(x,y+1) = y/(x+y) * B(x,y)")
print("-" * 70)
print(f"{'x':>8} {'y':>8} {'Left side':>18} {'Right side':>18} {'Error':>14}")

for x, y in test_values:
    lhs = beta(x, y + 1)
    rhs = (y / (x + y)) * beta(x, y)
    error = abs(lhs - rhs)
    print(f"{x:8.2f} {y:8.2f} {lhs:18.12f} {rhs:18.12f} {error:14.4e}")

print("\nBoth recurrence relations verified numerically.")


# ---------------------------------------------------------
# 5. Relation Between the Beta and Gamma Functions
#
# B(x,y) = Gamma(x) * Gamma(y) / Gamma(x+y)
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("RELATION BETWEEN THE BETA AND GAMMA FUNCTIONS")
print("=" * 70)

print("\nB(x,y) = Gamma(x) * Gamma(y) / Gamma(x+y)")
print("-" * 70)

print(f"{'x':>8} {'y':>8} {'Beta(x,y)':>18} "
      f"{'Gamma formula':>18} {'Error':>14}")

for x, y in test_values:
    beta_value = beta(x, y)

    gamma_value = (
        gamma(x) * gamma(y) / gamma(x + y)
    )

    error = abs(beta_value - gamma_value)

    print(
        f"{x:8.2f} {y:8.2f} "
        f"{beta_value:18.12f} "
        f"{gamma_value:18.12f} "
        f"{error:14.4e}"
    )

print("\nBeta-Gamma relation verified numerically.")

# ---------------------------------------------------------
# 6. Special Values of the Beta Function
#
# B(x,1) = 1/x
# B(1,y) = 1/y
#
# For positive integers m,n:
# B(m,n) = (m-1)! (n-1)! / (m+n-1)!
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("SPECIAL VALUES OF THE BETA FUNCTION")
print("=" * 70)

# Part 1: B(x,1) = 1/x

print("\n1. Verification of B(x,1) = 1/x")
print("-" * 70)
print(f"{'x':>8} {'Beta(x,1)':>18} {'1/x':>18} {'Error':>14}")

special_x_values = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]

for x in special_x_values:
    lhs = beta(x, 1)
    rhs = 1 / x
    error = abs(lhs - rhs)
    print(f"{x:8.2f} {lhs:18.12f} {rhs:18.12f} {error:14.4e}")


# Part 2: B(1,y) = 1/y

print("\n2. Verification of B(1,y) = 1/y")
print("-" * 70)
print(f"{'y':>8} {'Beta(1,y)':>18} {'1/y':>18} {'Error':>14}")

special_y_values = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]

for y in special_y_values:
    lhs = beta(1, y)
    rhs = 1 / y
    error = abs(lhs - rhs)
    print(f"{y:8.2f} {lhs:18.12f} {rhs:18.12f} {error:14.4e}")


# Part 3: Integer formula

print("\n3. Verification of the integer formula")
print("-" * 70)
print("B(m,n) = (m-1)! (n-1)! / (m+n-1)!")
print("-" * 70)

integer_pairs = [(1, 1), (1, 3), (2, 3), (3, 4), (4, 5), (5, 2)]

print(f"{'m':>6} {'n':>6} {'Beta(m,n)':>18} "
      f"{'Factorial formula':>20} {'Error':>14}")

for m, n in integer_pairs:
    lhs = beta(m, n)

    rhs = (
        math.factorial(m - 1) * math.factorial(n - 1)
        / math.factorial(m + n - 1)
    )

    error = abs(lhs - rhs)

    print(
        f"{m:6d} {n:6d} {lhs:18.12f} "
        f"{rhs:20.12f} {error:14.4e}"
    )

print("\nAll special values verified numerically.")

# ---------------------------------------------------------
# 7. Graphical Visualization of the Beta Function
# Part 1: B(x,y) for several fixed values of y
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("GRAPHICAL VISUALIZATION OF THE BETA FUNCTION")
print("=" * 70)

x_values = np.linspace(0.2, 5.0, 300)

fixed_y_values = [0.5, 1.0, 2.0, 3.0, 5.0]

plt.figure(figsize=(10, 6))

for y in fixed_y_values:
    beta_values = beta(x_values, y)

    plt.plot(
        x_values,
        beta_values,
        label=f"y = {y}"
    )

plt.xlabel("x", fontsize=12)
plt.ylabel("B(x,y)", fontsize=12)
plt.title("Beta Function for Several Fixed Values of y", fontsize=14)

plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("beta_function_curves.png", dpi=300)
plt.show()

print("Graph saved as beta_function_curves.png")

# ---------------------------------------------------------
# 7. Graphical Visualization of the Beta Function
# Part 2: Two-dimensional heatmap
# ---------------------------------------------------------

print("\nGenerating the Beta function heatmap...")

x_grid = np.linspace(0.2, 5.0, 200)
y_grid = np.linspace(0.2, 5.0, 200)

X, Y = np.meshgrid(x_grid, y_grid)

Z = beta(X, Y)

plt.figure(figsize=(9, 7))

contour = plt.contourf(X, Y, Z, levels=50, cmap="viridis")

plt.colorbar(contour, label="B(x,y)")

plt.xlabel("x", fontsize=12)
plt.ylabel("y", fontsize=12)
plt.title("Heatmap of the Beta Function", fontsize=14)

plt.tight_layout()

plt.savefig("beta_function_heatmap.png", dpi=300)
plt.show()

print("Heatmap saved as beta_function_heatmap.png")

# ---------------------------------------------------------
# 8. Three-Dimensional Surface Plot of the Beta Function
#
# z = B(x,y)
# ---------------------------------------------------------

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm

print("\nGenerating the 3D surface of the Beta function...")

x_3d = np.linspace(0.2, 5.0, 100)
y_3d = np.linspace(0.2, 5.0, 100)

X_3D, Y_3D = np.meshgrid(x_3d, y_3d)

Z_3D = beta(X_3D, Y_3D)

fig = plt.figure(figsize=(11, 8))
ax = fig.add_subplot(111, projection="3d")

surface = ax.plot_surface(
    X_3D,
    Y_3D,
    Z_3D,
    cmap="viridis",
    norm=LogNorm(),
    edgecolor="none",
    antialiased=True
)

ax.set_xlabel("x", labelpad=10)
ax.set_ylabel("y", labelpad=10)
ax.set_zlabel("B(x,y)", labelpad=10)

ax.set_title("Three-Dimensional Surface of the Beta Function")

fig.colorbar(
    surface,
    ax=ax,
    shrink=0.6,
    label="B(x,y) (logarithmic color scale)"
)

plt.tight_layout()

plt.savefig("beta_function_3d_surface.png", dpi=300)
plt.show()

print("3D surface saved as beta_function_3d_surface.png")

# ---------------------------------------------------------
# 9. Numerical Accuracy Comparison
#
# Method 1: SciPy beta function
# Method 2: Gamma-function quotient
# Method 3: Numerical integration
# ---------------------------------------------------------

print("\n" + "=" * 70)
print("NUMERICAL ACCURACY COMPARISON")
print("=" * 70)

print("\nComparison of three methods for evaluating B(x,y)")
print("-" * 70)

print(
    f"{'x':>6} {'y':>6} "
    f"{'SciPy beta':>18} "
    f"{'Gamma formula':>18} "
    f"{'Numerical integral':>20}"
)

for x, y in test_values:

    scipy_value = beta(x, y)

    gamma_value = (
        gamma(x) * gamma(y) / gamma(x + y)
    )

    integral_value, _ = quad(
        lambda t: t**(x - 1) * (1 - t)**(y - 1),
        0,
        1
    )

    print(
        f"{x:6.2f} {y:6.2f} "
        f"{scipy_value:18.12f} "
        f"{gamma_value:18.12f} "
        f"{integral_value:20.12f}"
    )

print("\nAbsolute errors relative to SciPy beta")
print("-" * 70)

print(
    f"{'x':>6} {'y':>6} "
    f"{'Gamma error':>18} "
    f"{'Integral error':>18}"
)

for x, y in test_values:

    scipy_value = beta(x, y)

    gamma_value = (
        gamma(x) * gamma(y) / gamma(x + y)
    )

    integral_value, _ = quad(
        lambda t: t**(x - 1) * (1 - t)**(y - 1),
        0,
        1
    )

    gamma_error = abs(scipy_value - gamma_value)
    integral_error = abs(scipy_value - integral_value)

    print(
        f"{x:6.2f} {y:6.2f} "
        f"{gamma_error:18.4e} "
        f"{integral_error:18.4e}"
    )

print("\nNumerical comparison completed.")
