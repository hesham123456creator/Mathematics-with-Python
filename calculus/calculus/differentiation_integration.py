
import sympy as sp

# Define the symbolic variable
x = sp.symbols('x', real=True)

# Define the function
f = sp.exp(x) + x**2 + 3*x

# First and second derivatives
first_derivative = sp.diff(f, x)
second_derivative = sp.diff(f, x, 2)

# Indefinite integral
integral = sp.integrate(f, x)

# Definite integral from 0 to 1
definite_integral = sp.integrate(f, (x, 0, 1))

# Verify the derivatives
first_check = sp.simplify(sp.diff(f, x) - first_derivative)
second_check = sp.simplify(sp.diff(f, x, 2) - second_derivative)

# Verify the indefinite integral
integral_check = sp.simplify(sp.diff(integral, x) - f)

# Numerical evaluation at x = 1
numerical_value = sp.N(f.subs(x, 1))

# Display results
print("Function:", f)
print("First derivative:", first_derivative)
print("Second derivative:", second_derivative)
print("Indefinite integral:", integral)
print("Definite integral from 0 to 1:", definite_integral)

print("\nVerification:")
print("First derivative verified:", first_check == 0)
print("Second derivative verified:", second_check == 0)
print("Integral verified:", integral_check == 0)

print("\nNumerical evaluation:")
print("f(1) =", numerical_value)
