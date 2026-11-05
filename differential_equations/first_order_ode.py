import sympy as sp

# Define the independent variable and unknown function
x = sp.symbols('x')
y = sp.Function('y')

# Define the differential equation:
# y' + 2y = exp(x)
ode = sp.Eq(sp.diff(y(x), x) + 2*y(x), sp.exp(x))

# Solve the differential equation
solution = sp.dsolve(ode, y(x))

print("Differential equation:")
sp.pprint(ode)

print("\nGeneral solution:")
sp.pprint(solution)
