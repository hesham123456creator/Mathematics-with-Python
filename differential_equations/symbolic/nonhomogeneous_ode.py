
"""
Solve and visualize a nonhomogeneous second-order ODE using SymPy and SciPy.

The script derives the symbolic solution, verifies initial conditions,
and compares the exact solution with a numerical RK45 approximation.
"""

from pathlib import Path

# Pylint does not fully recognize SymPy symbolic functions.
# pylint: disable=not-callable

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Directory for generated plots
IMAGE_DIR = Path(__file__).resolve().parents[2] / "images"
IMAGE_DIR.mkdir(parents=True, exist_ok=True)


# Define symbolic variables and unknown function
x = sp.symbols("x", real=True)
C1, C2 = sp.symbols("C1 C2")
y = sp.Function("y")


# Nonhomogeneous second-order ODE:
# y'' - 3y' + 2y = exp(x)
ode = sp.Eq(
    sp.diff(y(x), x, 2)
    - 3 * sp.diff(y(x), x)
    + 2 * y(x),
    sp.exp(x)
)


# Solve the differential equation
solution = sp.dsolve(ode, y(x))
general_solution = solution.rhs


# Verify the general solution
verification = sp.simplify(
    sp.diff(general_solution, x, 2)
    - 3 * sp.diff(general_solution, x)
    + 2 * general_solution
    - sp.exp(x)
)


# Display the equation and general solution
print("NONHOMOGENEOUS SECOND-ORDER ODE")
print("=" * 60)

print("\nDifferential equation:")
sp.pprint(ode)

print("\nGeneral solution:")
sp.pprint(solution)

print("\nSymbolic verification:")
sp.pprint(verification)

if verification == 0:
    print("Verification successful: the solution satisfies the ODE.")
else:
    print("Verification failed.")


# Initial-value problem:
# y(0) = 1
# y'(0) = 0

initial_conditions = [
    sp.Eq(general_solution.subs(x, 0), 1),
    sp.Eq(sp.diff(general_solution, x).subs(x, 0), 0)
]

constants = sp.solve(initial_conditions, (C1, C2))

particular_solution = sp.simplify(
    general_solution.subs(constants)
)


# Display initial-value solution
print("\nINITIAL-VALUE PROBLEM")
print("=" * 60)
print("Initial conditions: y(0) = 1, y'(0) = 0")

print("\nIntegration constants:")
sp.pprint(constants)

print("\nParticular solution:")
sp.pprint(particular_solution)


# Verify initial conditions
y_initial = sp.simplify(
    particular_solution.subs(x, 0)
)

dy_initial = sp.simplify(
    sp.diff(particular_solution, x).subs(x, 0)
)

y_residual = sp.simplify(y_initial - 1)
dy_residual = sp.simplify(dy_initial)

print("\nVerification of initial conditions:")
print("y(0) =", y_initial)
print("y'(0) =", dy_initial)

print("\nInitial-condition residuals:")
print("y(0) - 1 =", y_residual)
print("y'(0) - 0 =", dy_residual)

if y_residual == 0 and dy_residual == 0:
    print("Initial conditions verified successfully.")
else:
    print("Initial-condition verification failed.")


# Verify particular solution
ivp_verification = sp.simplify(
    sp.diff(particular_solution, x, 2)
    - 3 * sp.diff(particular_solution, x)
    + 2 * particular_solution
    - sp.exp(x)
)

print("\nVerification of particular solution:")
sp.pprint(ivp_verification)

if ivp_verification == 0:
    print("Particular solution satisfies the ODE.")
else:
    print("Particular solution verification failed.")


# Convert symbolic expressions into numerical functions
x_values = np.linspace(0, 3, 300)

general_function = sp.lambdify(
    (x, C1, C2),
    general_solution,
    modules="numpy"
)

particular_function = sp.lambdify(
    x,
    particular_solution,
    modules="numpy"
)

# Numerical solution using SciPy's RK45 method


def ode_system(x_value, state):
    """Define the first-order system equivalent to the second-order ODE."""
    y_value, dy_value = state

    return [
        dy_value,
        3 * dy_value - 2 * y_value + np.exp(x_value)
    ]


numerical_solution = solve_ivp(
    ode_system,
    (0, 3),
    [1, 0],
    method="RK45",
    t_eval=x_values,
    rtol=1e-9,
    atol=1e-11
)

if not numerical_solution.success:
    raise RuntimeError(numerical_solution.message)


# Compare numerical and analytical solutions

exact_values = particular_function(x_values)
numerical_values = numerical_solution.y[0]

absolute_errors = np.abs(
    numerical_values - exact_values
)

maximum_error = np.max(absolute_errors)

print("\nNUMERICAL VERIFICATION")
print("=" * 60)
print("Numerical method: RK45")
print("Number of evaluation points:", len(x_values))
print(f"Maximum absolute error: {maximum_error:.6e}")


# Plot several general solutions
plt.figure(figsize=(9, 5))

for c1, c2 in [(0, 0), (1, 0), (0, 1), (1, -1), (-1, 1)]:
    plt.plot(
        x_values,
        general_function(x_values, c1, c2),
        label=f"C1={c1}, C2={c2}"
    )

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("General Solutions of a Nonhomogeneous Second-Order ODE")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
    IMAGE_DIR / "nonhomogeneous_ode_general.png",
    dpi=300
)

plt.show()


# Plot the particular solution
plt.figure(figsize=(9, 5))

plt.plot(
    x_values,
    particular_function(x_values),
    label="Particular solution",
    linewidth=2
)

plt.scatter(
    [0],
    [1],
    label="Initial point (0, 1)"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonhomogeneous ODE: Initial-Value Solution")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
    IMAGE_DIR / "nonhomogeneous_ode.png",
    dpi=300
)


# Compare analytical and numerical solutions

plt.figure(figsize=(9, 5))

plt.plot(
    x_values,
    exact_values,
    label="Analytical solution",
    linewidth=2
)

plt.plot(
    numerical_solution.t,
    numerical_values,
    "--",
    label="Numerical solution (RK45)",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Analytical vs Numerical Solution")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
    IMAGE_DIR / "nonhomogeneous_ode_comparison.png",
    dpi=300
)


plt.show()
