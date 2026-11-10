
# Pylint does not fully recognize SymPy symbolic functions.
# pylint: disable=not-callable

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


# Define symbolic variables and unknown function
x = sp.symbols("x", real=True)
C1 = sp.symbols("C1")
y = sp.Function("y")

# Define the first-order linear ODE:
# y' + 2y = exp(x)
ode = sp.Eq(
    sp.diff(y(x), x) + 2*y(x),
    sp.exp(x)
)


# Solve the differential equation
solution = sp.dsolve(ode, y(x))

# Extract the right-hand side
general_solution = solution.rhs


# Verify the symbolic solution
verification = sp.simplify(
    sp.diff(general_solution, x)
    + 2*general_solution
    - sp.exp(x)
)


# Display the differential equation and solution
print("FIRST-ORDER LINEAR DIFFERENTIAL EQUATION")
print("=" * 55)

print("\nDifferential equation:")
sp.pprint(ode)

print("\nGeneral solution:")
sp.pprint(solution)

print("\nSymbolic verification:")
print("Substitution into the ODE gives:")
sp.pprint(verification)

if verification == 0:
    print("Verification successful: the solution satisfies the ODE.")
else:
    print("Verification failed.")


# Initial-value problem:
# y(0) = 1
initial_condition = {x: 0}

constant = sp.solve(
    general_solution.subs(x, 0) - 1,
    C1
)[0]

particular_solution = sp.simplify(
    general_solution.subs(C1, constant)
)

print("\nINITIAL-VALUE PROBLEM")
print("=" * 55)
print("Initial condition: y(0) = 1")

print("\nConstant of integration:")
sp.pprint(constant)

print("\nParticular solution:")
sp.pprint(particular_solution)

# Verify initial condition
initial_value = particular_solution.subs(x, 0)

print("\nVerification of initial condition:")
sp.pprint(initial_value)


# Numerical evaluation for plotting
x_values = np.linspace(-2, 2, 300)

general_function = sp.lambdify(
    (x, C1),
    general_solution,
    modules="numpy"
)

particular_function = sp.lambdify(
    x,
    particular_solution,
    modules="numpy"
)

# Plot several members of the general solution
plt.figure(figsize=(9, 5))

for c in [-2, -1, 0, 1, 2]:
    plt.plot(
        x_values,
        general_function(x_values, c),
        label=f"C = {c}"
    )

# Plot the particular solution
plt.plot(
    x_values,
    particular_function(x_values),
    "k--",
    linewidth=2,
    label="Solution with y(0) = 1"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Symbolic Solution of a First-Order Linear ODE")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("first_order_ode.png", dpi=300)
plt.show()
