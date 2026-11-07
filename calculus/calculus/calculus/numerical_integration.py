
import numpy as np
import sympy as sp
from scipy.integrate import quad, trapezoid, simpson

# Define the integrand
def f(x):
    return np.exp(-x**2)

# Exact integral
x = sp.symbols('x')
exact_value = float(
    sp.sqrt(sp.pi) / 2 * sp.erf(1)
)

# Adaptive quadrature
quad_result, quad_error = quad(f, 0, 1)

# Numerical integration using equally spaced points
n = 100
x_values = np.linspace(0, 1, n + 1)
y_values = f(x_values)

trapezoid_result = trapezoid(y_values, x_values)
simpson_result = simpson(y_values, x=x_values)

# Display results
print("Numerical Integration of exp(-x^2) from 0 to 1")
print("-" * 55)
print(f"Exact value:       {exact_value:.12f}")
print(f"Quad result:       {quad_result:.12f}")
print(f"Trapezoid result:  {trapezoid_result:.12f}")
print(f"Simpson result:    {simpson_result:.12f}")
print("-" * 55)

print("Absolute errors:")
print(f"Quad:      {abs(quad_result - exact_value):.3e}")
print(f"Trapezoid: {abs(trapezoid_result - exact_value):.3e}")
print(f"Simpson:   {abs(simpson_result - exact_value):.3e}")
print(f"Quad estimated error: {quad_error:.3e}")
