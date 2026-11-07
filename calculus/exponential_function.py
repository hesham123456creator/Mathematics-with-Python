
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# 1. Define the domain
# --------------------------------------------------

x = np.linspace(-3, 3, 400)

# --------------------------------------------------
# 2. Define the exponential function and derivative
# --------------------------------------------------

f = np.exp(x)
f_prime = np.exp(x)

# Tangent line at x = 0
tangent = 1 + x

# --------------------------------------------------
# 3. Numerical verification
# --------------------------------------------------

numerical_derivative = np.gradient(f, x)

derivative_error = np.abs(
    numerical_derivative - f_prime
)

print("=" * 55)
print("EXPONENTIAL FUNCTION: NUMERICAL ANALYSIS")
print("=" * 55)

print("\nMathematical property:")
print("f(x) = exp(x)")
print("f'(x) = exp(x)")
print("Therefore, f'(x) = f(x).")

print("\nNumerical verification:")
print(
    "Maximum derivative approximation error:",
    f"{np.max(derivative_error):.8e}"
)

# --------------------------------------------------
# 4. Plot the exponential function and tangent
# --------------------------------------------------

plt.figure(figsize=(10, 6))

plt.plot(
    x, f,
    label=r"$f(x)=e^x$",
    linewidth=2
)

plt.plot(
    x, tangent,
    "--",
    label=r"Tangent at $x=0$: $y=1+x$"
)

plt.scatter(
    [0], [1],
    label="Point of tangency (0, 1)"
)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Exponential Function and Tangent Line")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig(
    "exponential_function.png",
    dpi=300
)

plt.show()

# --------------------------------------------------
# 5. Compare function and derivative
# --------------------------------------------------

plt.figure(figsize=(10, 6))

plt.plot(
    x, f,
    label=r"$f(x)=e^x$",
    linewidth=2
)

plt.plot(
    x, f_prime,
    "--",
    label=r"$f'(x)=e^x$",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("Value")
plt.title("Exponential Function and Its Derivative")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig(
    "exponential_derivative.png",
    dpi=300
)

plt.show()

# --------------------------------------------------
# 6. Logarithmic-scale visualization
# --------------------------------------------------

plt.figure(figsize=(10, 6))

plt.semilogy(
    x, f,
    label=r"$f(x)=e^x$",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("f(x), logarithmic scale")
plt.title("Exponential Growth on a Logarithmic Scale")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()

plt.savefig(
    "exponential_function_log_scale.png",
    dpi=300
)

plt.show()

print("\nGenerated figures:")
print("- exponential_function.png")
print("- exponential_derivative.png")
print("- exponential_function_log_scale.png")

print("\nAnalysis completed successfully.")
