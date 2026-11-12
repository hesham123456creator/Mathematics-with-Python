
"""
Solve a second-order boundary value problem using SciPy.

The script solves y'' + y = 0 subject to two boundary conditions,
compares the numerical solution with the exact solution y(x) = sin(x),
and verifies the numerical accuracy.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp


# --------------------------------------------------
# Boundary Value Problem
#
# y'' + y = 0
# y(0) = 0
# y(pi/2) = 1
#
# Exact solution: y(x) = sin(x)
# --------------------------------------------------

# Save generated plots in the repository images directory.
IMAGE_DIR = Path(__file__).resolve().parents[2] / "images"
IMAGE_DIR.mkdir(parents=True, exist_ok=True)


# 1. Convert the second-order equation into a
#    first-order system:
#    y1' = y2
#    y2' = -y1

def ode(x, y):
    """Return the first-order system."""
    return np.vstack((y[1], -y[0]))


# 2. Boundary conditions

def boundary_conditions(ya, yb):
    """Return the residuals of the boundary conditions."""
    return np.array([ya[0], yb[0] - 1.0])


# 3. Initial mesh and guess

a = 0.0
b = np.pi / 2

x = np.linspace(a, b, 20)
y_guess = np.zeros((2, x.size))


# 4. Solve using SciPy solve_bvp

solution = solve_bvp(
    ode,
    boundary_conditions,
    x,
    y_guess,
    tol=1e-8,
    max_nodes=1000,
)

if not solution.success:
    raise RuntimeError(f"BVP solver failed: {solution.message}")


# 5. Evaluate the numerical and exact solutions

x_plot = np.linspace(a, b, 200)

y_numerical = solution.sol(x_plot)[0]
y_exact = np.sin(x_plot)


# 6. Error analysis

absolute_errors = np.abs(y_numerical - y_exact)
maximum_error = np.max(absolute_errors)

left_residual = abs(y_numerical[0])
right_residual = abs(y_numerical[-1] - 1.0)

boundary_tolerance = 1e-8

if max(left_residual, right_residual) > boundary_tolerance:
    raise RuntimeError("Numerical boundary conditions were not satisfied.")


# 7. Display numerical results

print("Boundary Value Problem Using SciPy solve_bvp")
print("=" * 60)

print("Equation: y'' + y = 0")
print("Boundary conditions: y(0) = 0, y(pi/2) = 1")

print("\nSolver Information")
print("Solver successful:", solution.success)
print("Solver message:", solution.message)
print("Number of mesh points:", solution.x.size)

print("\nAccuracy Analysis")
print(f"Maximum absolute error: {maximum_error:.10e}")

print("\nBoundary Condition Verification")
print(f"y(0): {y_numerical[0]:.10e}")
print(f"y(pi/2): {y_numerical[-1]:.10e}")
print(f"Left boundary residual: {left_residual:.10e}")
print(f"Right boundary residual: {right_residual:.10e}")

print("\nSelected Numerical Values")
print("-" * 45)
print(f"{'x':>12} {'Numerical':>15} {'Exact':>15}")
for x_value in np.linspace(a, b, 5):
    numerical_value = solution.sol(x_value)[0]
    exact_value = np.sin(x_value)
    print(f"{x_value:12.6f} {numerical_value:15.10f} {exact_value:15.10f}")


# 8. Plot numerical and exact solutions

plt.figure(figsize=(9, 6))

plt.plot(
    x_plot,
    y_numerical,
    label="Numerical solution",
    linewidth=2.5,
)

plt.plot(
    x_plot,
    y_exact,
    "--",
    label="Exact solution",
    linewidth=2,
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Boundary Value Problem Using SciPy solve_bvp")

plt.legend()
plt.grid(True)
plt.tight_layout()

plot_path = IMAGE_DIR / "solve_bvp_example.png"
plt.savefig(plot_path, dpi=300)

print(f"\nPlot saved to: {plot_path}")

plt.show()
plt.close()
