
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.optimize import brentq

# --------------------------------------------------
# Nonlinear Bratu Boundary Value Problem
#
# y'' + lambda * exp(y) = 0
# y(0) = 0, y(1) = 0
# --------------------------------------------------

lam = 1.0

# --------------------------------------------------
# 1. First-order system
# --------------------------------------------------

def ode(x, y):
    return np.vstack((
        y[1],
        -lam * np.exp(y[0])
    ))


# Boundary conditions
def bc(ya, yb):
    return np.array([
        ya[0],
        yb[0]
    ])


# --------------------------------------------------
# 2. Exact analytical solution
# --------------------------------------------------

# The Bratu equation has the analytical solution:
#
# y(x) = 2 log(cosh(a/2) / cosh(a*(x-1/2)))
#
# where lambda = 2*a^2 / cosh(a/2)^2.
#
# We select the smaller positive root.

def parameter_equation(a):
    return 2 * a**2 / np.cosh(a / 2)**2 - lam


a = brentq(parameter_equation, 1e-8, 2.4)


def exact_solution(x):
    return 2 * np.log(
        np.cosh(a / 2) /
        np.cosh(a * (x - 0.5))
    )


# --------------------------------------------------
# 3. Initial mesh and initial guess
# --------------------------------------------------

x = np.linspace(0.0, 1.0, 30)

y_guess = np.zeros((2, x.size))
y_guess[0] = 0.1 * np.sin(np.pi * x)


# --------------------------------------------------
# 4. Solve using SciPy solve_bvp
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
# 5. Evaluate numerical and exact solutions
# --------------------------------------------------

x_plot = np.linspace(0.0, 1.0, 300)

y_numerical = solution.sol(x_plot)[0]
y_exact = exact_solution(x_plot)

# Maximum absolute error
maximum_error = np.max(
    np.abs(y_numerical - y_exact)
)

# Boundary residuals
left_residual = abs(y_numerical[0])
right_residual = abs(y_numerical[-1])

# ODE residual at plotting points
solution_derivative = solution.sol(x_plot, 1)

ode_residual = np.max(
    np.abs(
        solution_derivative[1]
        + lam * np.exp(y_numerical)
    )
)

scipy_rms_residual = np.max(solution.rms_residuals)

# --------------------------------------------------
# 6. Display numerical results
# --------------------------------------------------

print("Nonlinear Bratu Boundary Value Problem")
print("=" * 60)

print(f"Parameter lambda: {lam}")
print(f"Analytical parameter a: {a:.10f}")

print("\nSolver Information")
print("Solver successful:", solution.success)
print("Solver message:", solution.message)
print("Number of mesh points:", solution.x.size)

print("\nAccuracy Analysis")
print(f"Maximum solution value: {np.max(y_numerical):.10f}")
print(f"Maximum absolute error: {maximum_error:.10e}")
print(f"ODE residual (plot points): {ode_residual:.10e}")
print(f"Maximum SciPy RMS residual: {scipy_rms_residual:.10e}")

print("\nBoundary Conditions")
print(f"y(0): {y_numerical[0]:.10e}")
print(f"y(1): {y_numerical[-1]:.10e}")
print(f"Left boundary residual: {left_residual:.10e}")
print(f"Right boundary residual: {right_residual:.10e}")


# --------------------------------------------------
# 7. Plot numerical and exact solutions
# --------------------------------------------------

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
    label="Exact analytical solution",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonlinear Bratu Boundary Value Problem")

plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
"images/nonlinear_bvp_nontrivial.png",
    dpi=300
)

plt.show()
plt.close()

