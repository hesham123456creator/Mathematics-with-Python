
import sympy as sp

# --------------------------------------------------
# 1. Define the symbolic variable and function
# --------------------------------------------------

x = sp.symbols("x", real=True)

f = sp.exp(x) + x**2 + 3*x

print("=" * 60)
print("SYMBOLIC DIFFERENTIATION AND INTEGRATION")
print("=" * 60)

print("\nOriginal function:")
sp.pprint(f)

# --------------------------------------------------
# 2. First and second derivatives
# --------------------------------------------------

first_derivative = sp.diff(f, x)
second_derivative = sp.diff(f, x, 2)

print("\nFirst derivative:")
sp.pprint(first_derivative)

print("\nSecond derivative:")
sp.pprint(second_derivative)

# --------------------------------------------------
# 3. Symbolic verification of derivatives
# --------------------------------------------------

verification_first = sp.simplify(
    sp.diff(f, x) - first_derivative
)

verification_second = sp.simplify(
    sp.diff(f, x, 2) - second_derivative
)

print("\nDerivative verification:")
print("First derivative verified:", verification_first == 0)
print("Second derivative verified:", verification_second == 0)

# --------------------------------------------------
# 4. Indefinite integration
# --------------------------------------------------

indefinite_integral = sp.integrate(f, x)

print("\nIndefinite integral:")
sp.pprint(indefinite_integral)

# Verify by differentiation
integral_verification = sp.simplify(
    sp.diff(indefinite_integral, x) - f
)

print(
    "Integration verified:",
    integral_verification == 0
)

# --------------------------------------------------
# 5. Definite integration
# --------------------------------------------------

a = 0
b = 1

definite_integral = sp.integrate(f, (x, a, b))

print(f"\nDefinite integral from {a} to {b}:")
sp.pprint(definite_integral)

print(
    "Numerical value:",
    sp.N(definite_integral, 10)
)

# --------------------------------------------------
# 6. Numerical evaluation
# --------------------------------------------------

x_value = 1

print(f"\nFunction evaluation at x = {x_value}:")
print("f(x) =", sp.N(f.subs(x, x_value), 10))

print("\nFirst derivative at x = 1:")
print(
    sp.N(first_derivative.subs(x, x_value), 10)
)

print("\nSecond derivative at x = 1:")
print(
    sp.N(second_derivative.subs(x, x_value), 10)
)

print("\n" + "=" * 60)
print("COMPUTATION COMPLETED SUCCESSFULLY")
print("=" * 60)
