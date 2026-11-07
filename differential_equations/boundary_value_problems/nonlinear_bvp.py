
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# --------------------------------------------------
# Nonlinear Boundary Value Problem
#
# y'' + y^2 = 0
# y(0) = 0, y(1) = 0
#
# The zero function is an exact solution.
# --------------------------------------------------

# Convert the second-order equation into a
# first-order system:
# y1' = y2
# y2' = -y1^2

def ode(x, y):
    return np.vstack((
        y[1],
        -y[0]**2
    ))


# Boundary conditions
def bc(ya, yb):
    return np.array([
        ya[0],
        yb[0]
    ])


# --------------------------------------------------
# 1. Initial mesh and initial guess
# --------------------------------------------------

x = np.linspace(0.0, 1.0, 30)

# Initial guess corresponding to the trivial solution
y_guess = np.zeros((2, x.size))


# --------------------------------------------------
# 2. Solve the nonlinear BVP
# --------------------------------------------------

solution = solve_bvp(
    ode,
    bc,
    x,
    y_guess,
    tol=1e-8,
    max_nodes=1000
)

if not solution.success:
    raise RuntimeError(
        "BVP solver failed: " + solution.message
    )


# --------------------------------------------------
# 3. Evaluate the numerical solution
# --------------------------------------------------

x_plot = np.linspace(0.0, 1.0, 200)
y_plot = solution.sol(x_plot)[0]

# Exact trivial solution
y_exact = np.zeros_like(x_plot)

# Maximum absolute error
error = np.max(np.abs(y_plot - y_exact))

# Numerical boundary residuals
left_residual = abs(y_plot[0])
right_residual = abs(y_plot[-1])


# --------------------------------------------------
# 4. Display numerical results
# --------------------------------------------------

print("Nonlinear Boundary Value Problem")
print("=" * 55)

print("Equation: y'' + y^2 = 0")
print("Boundary conditions: y(0) = 0, y(1) = 0")

print("\nSolver successful:", solution.success)
print("Solver message:", solution.message)
print("Number of mesh points:", solution.x.size)

print("\nNumerical Results")
print(f"Maximum solution value: {np.max(np.abs(y_plot)):.10e}")
print(f"Maximum absolute error: {error:.10e}")
print(f"Left boundary residual: {left_residual:.10e}")
print(f"Right boundary residual: {right_residual:.10e}")


# --------------------------------------------------
# 5. Plot numerical and exact solutions
# --------------------------------------------------

plt.figure(figsize=(9, 6))

plt.plot(
    x_plot,
    y_plot,
    label="Numerical solution",
    linewidth=2
)

plt.plot(
    x_plot,
    y_exact,
    "--",
    label="Exact trivial solution",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonlinear Boundary Value Problem")

plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("nonlinear_bvp.png", dpi=300)

plt.show()
plt.close()

