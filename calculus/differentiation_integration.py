
import sympy as sp

# Define the symbolic variable
x = sp.symbols('x')

# Define the function
f = sp.exp(x) + x**2 + 3*x

# First derivative
first_derivative = sp.diff(f, x)

# Second derivative
second_derivative = sp.diff(f, x, 2)

# Indefinite integral
integral = sp.integrate(f, x)

# Display results
print("Function:", f)
print("First derivative:", first_derivative)
print("Second derivative:", second_derivative)
print("Indefinite integral:", integral)
