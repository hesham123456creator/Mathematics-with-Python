
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Boundary-value problem:
# y'' = -pi^2 sin(pi*x)
# y(0) = 0, y(1) = 0
# Exact solution: y(x) = sin(pi*x)

def system(x, y):
    return np.vstack((
        y[1],
        -np.pi**2 * np.sin(np.pi * x)
    ))


def boundary_conditions(ya, yb):
    return np.array([ya[0], yb[0]])


# Initial mesh and guess
x = np.linspace(0.0, 1.0, 20)
y_guess = np.zeros((2, x.size))

# Solve the BVP
solution = solve_bvp(
    system,
    boundary_conditions,
    x,
    y_guess,
    tol=1e-8,
    max_nodes=1000
)

if not solution.success:
    raise RuntimeError(solution.message)

# Evaluate numerical solution
x_plot = np.linspace(0.0, 1.0, 300)
y_numerical = solution.sol(x_plot)[0]
y_derivative = solution.sol(x_plot)[1]

# Exact solution and derivative
y_exact = np.sin(np.pi * x_plot)
derivative_exact = np.pi * np.cos(np.pi * x_plot)

# Error analysis
absolute_error = np.abs(y_numerical - y_exact)
derivative_error = np.abs(y_derivative - derivative_exact)

# Boundary-condition residuals
boundary_residuals = boundary_conditions(
    solution.sol(0.0),
    solution.sol(1.0)
)

# Differential-equation residual
ode_residual = (
    solution.sol(x_plot, 2)[0]
    + np.pi**2 * np.sin(np.pi * x_plot)
)

print("SciPy Boundary Value Problem Solver")
print("-" * 50)
print("Solver successful:", solution.success)
print("Solver message:", solution.message)
print("Number of mesh nodes:", solution.x.size)
print("Maximum solution error:", f"{np.max(absolute_error):.6e}")
print("Maximum derivative error:", f"{np.max(derivative_error):.6e}")
print("Maximum ODE residual:", f"{np.max(np.abs(ode_residual)):.6e}")
print("Boundary residuals:", boundary_residuals)

# Plot numerical and exact solutions
plt.figure(figsize=(9, 5))
plt.plot(x_plot, y_numerical, label="Numerical solution")
plt.plot(x_plot, y_exact, "--", label="Exact solution")
plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Boundary Value Problem Using solve_bvp")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("solve_bvp_example.png", dpi=300)
plt.show()

# Plot absolute error
plt.figure(figsize=(9, 5))
plt.semilogy(x_plot, absolute_error, label="Solution error")
plt.semilogy(x_plot, derivative_error, label="Derivative error")
plt.xlabel("x")
plt.ylabel("Absolute error")
plt.title("Numerical Error Analysis")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("solve_bvp_error.png", dpi=300)
plt.show()

# Plot differential-equation residual
plt.figure(figsize=(9, 5))
plt.plot(x_plot, ode_residual, label="ODE residual")
plt.xlabel("x")
plt.ylabel("Residual")
plt.title("Differential Equation Residual")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("solve_bvp_residual.png", dpi=300)
plt.show()
