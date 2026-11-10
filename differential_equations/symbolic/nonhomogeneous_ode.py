# Pylint does not fully recognize SymPy symbolic functions.
# pylint: disable=not-callable

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


# Define symbolic variables and unknown function
x = sp.symbols("x", real=True)
C1, C2 = sp.symbols("C1 C2")
y = sp.Function("y")


# Nonhomogeneous second-order ODE:
# y'' - 3y' + 2y = exp(x)
ode = sp.Eq(
    sp.diff(y(x), x, 2)
    - 3*sp.diff(y(x), x)
    + 2*y(x),
    sp.exp(x)
)


# Solve the differential equation
solution = sp.dsolve(ode, y(x))
general_solution = solution.rhs


# Verify the symbolic solution
verification = sp.simplify(
    sp.diff(general_solution, x, 2)
    - 3*sp.diff(general_solution, x)
    + 2*general_solution
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
y_initial = particular_solution.subs(x, 0)
dy_initial = sp.diff(particular_solution, x).subs(x, 0)

print("\nVerification of initial conditions:")
print("y(0) =", y_initial)
print("y'(0) =", dy_initial)


# Verify particular solution
ivp_verification = sp.simplify(
    sp.diff(particular_solution, x, 2)
    - 3*sp.diff(particular_solution, x)
    + 2*particular_solution
    - sp.exp(x)
)

print("\nVerification of particular solution:")
sp.pprint(ivp_verification)


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
plt.savefig("nonhomogeneous_ode_general.png", dpi=300)
plt.show()


# Plot the particular solution
plt.figure(figsize=(9, 5))

plt.plot(
    x_values,
    particular_function(x_values),
    label="Particular solution",
    linewidth=2
)

plt.scatter([0], [1], label="Initial point (0, 1)")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonhomogeneous ODE: Initial-Value Solution")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("nonhomogeneous_ode.png", dpi=300)
plt.show()
