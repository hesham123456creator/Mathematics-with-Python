
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

# --------------------------------------------------
# Boundary Value Problem
# y''(x) = -pi^2 sin(pi*x)
# y(0) = 0, y(1) = 0
# Exact solution: y(x) = sin(pi*x)
# --------------------------------------------------

# Define the interval and discretization
a = 0.0
b = 1.0
N = 40

h = (b - a) / N
x = np.linspace(a, b, N + 1)

# --------------------------------------------------
# 1. Construct the finite difference system
# --------------------------------------------------

diagonal = -2.0 * np.ones(N - 1)
off_diagonal = np.ones(N - 2)

# Banded representation of the tridiagonal matrix
ab = np.zeros((3, N - 1))
ab[0, 1:] = off_diagonal
ab[1, :] = diagonal
ab[2, :-1] = off_diagonal

# Right-hand side
rhs = -h**2 * np.pi**2 * np.sin(np.pi * x[1:N])

# --------------------------------------------------
# 2. Solve the linear system
# --------------------------------------------------

y_interior = solve_banded((1, 1), ab, rhs)

# Include boundary values
y_numerical = np.zeros(N + 1)
y_numerical[1:N] = y_interior

# --------------------------------------------------
# 3. Exact solution
# --------------------------------------------------

y_exact = np.sin(np.pi * x)

# --------------------------------------------------
# 4. Error analysis
# --------------------------------------------------

absolute_error = np.abs(y_numerical - y_exact)

maximum_error = np.max(absolute_error)
relative_error = maximum_error / np.max(np.abs(y_exact))

print("Finite Difference Method for a BVP")
print("=" * 55)

print(f"Number of subintervals: {N}")
print(f"Step size: {h:.8f}")
print(f"Maximum absolute error: {maximum_error:.10e}")
print(f"Maximum relative error: {relative_error:.10e}")

# --------------------------------------------------
# 5. Plot numerical and exact solutions
# --------------------------------------------------

plt.figure(figsize=(9, 6))

plt.plot(
    x,
    y_numerical,
    "o-",
    label="Finite Difference Solution",
    markersize=4
)

plt.plot(
    x,
    y_exact,
    "--",
    label="Exact Solution",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Finite Difference Method for a Boundary Value Problem")

plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("finite_difference_method.png", dpi=300)

plt.show()
plt.close()

