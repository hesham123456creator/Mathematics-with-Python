import sympy as sp

# Define the independent variable and unknown function
x = sp.symbols('x')
y = sp.Function('y')

# Define the second-order differential equation:
# y'' - 3y' + 2y = 0
ode = sp.Eq(
    sp.diff(y(x), x, 2) - 3*sp.diff(y(x), x) + 2*y(x),
    0
)

# Solve the differential equation
solution = sp.dsolve(ode, y(x))

print("Differential equation:")
sp.pprint(ode)

print("\nGeneral solution:")
sp.pprint(solution)
