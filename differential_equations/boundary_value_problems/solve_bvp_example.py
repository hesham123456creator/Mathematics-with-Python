
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

# 1. Convert the second-order equation into
#    a first-order system:
#    y1' = y2
#    y2' = -y1

def ode(x, y):
    return np.vstack((
        y[1],
        -y[0]
    ))


# 2. Boundary conditions

def boundary_conditions(ya, yb):
    return np.array([
        ya[0],
        yb[0] - 1.0
    ])


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
    max_nodes=1000
)

if not solution.success:
    raise RuntimeError(
        "BVP solver failed: " + solution.message
    )


# 5. Evaluate the numerical solution

x_plot = np.linspace(a, b, 200)

y_numerical = solution.sol(x_plot)[0]

# Exact solution
y_exact = np.sin(x_plot)


# 6. Error analysis

absolute_errors = np.abs(y_numerical - y_exact)

maximum_error = np.max(absolute_errors)

# Boundary residuals
left_residual = abs(y_numerical[0])
right_residual = abs(y_numerical[-1] - 1.0)


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


# 8. Plot numerical and exact solutions

plt.figure(figsize=(9, 6))

plt.plot(
    x_plot,
    y_numerical,
    label="Numerical solution",
    linewidth=2.5
)

plt.plot(
    x_plot,
    y_exact,
    "--",
    label="Exact solution",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Boundary Value Problem Using SciPy solve_bvp")

plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("solve_bvp_example.png", dpi=300)

plt.show()
plt.close()

