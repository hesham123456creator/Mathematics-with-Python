
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

# Boundary-value problem:
# y'' = -pi^2 sin(pi*x)
# y(0) = 0, y(1) = 0
# Exact solution: y(x) = sin(pi*x)

def solve_finite_difference(n):
    """Solve the BVP using the finite difference method."""

    h = 1.0 / n
    x = np.linspace(0.0, 1.0, n + 1)
    x_interior = x[1:-1]

    # Construct the tridiagonal coefficient matrix
    ab = np.zeros((3, n - 1))
    ab[0, 1:] = 1.0
    ab[1, :] = -2.0
    ab[2, :-1] = 1.0
    ab /= h**2

    # Right-hand side
    b = -np.pi**2 * np.sin(np.pi * x_interior)

    # Solve the banded linear system
    y_interior = solve_banded((1, 1), ab, b)

    # Include boundary conditions
    y = np.zeros(n + 1)
    y[1:-1] = y_interior

    return x, y


# Main computation
n = 20
x, y_numerical = solve_finite_difference(n)

# Exact solution
y_exact = np.sin(np.pi * x)

# Error analysis
absolute_error = np.abs(y_numerical - y_exact)
max_error = np.max(absolute_error)
mean_error = np.mean(absolute_error)
rms_error = np.sqrt(np.mean(absolute_error**2))

print("Finite Difference Method")
print("-" * 45)
print("Number of subintervals:", n)
print("Step size:", 1.0 / n)
print(f"Maximum absolute error: {max_error:.6e}")
print(f"Mean absolute error:    {mean_error:.6e}")
print(f"RMS error:              {rms_error:.6e}")

# Convergence study
grid_sizes = [10, 20, 40, 80, 160, 320]
errors = []

print("\nConvergence Study")
print("-" * 45)
print(f"{'n':>8} {'Maximum error':>20} {'Order':>12}")

for i, grid_n in enumerate(grid_sizes):
    grid_x, grid_y = solve_finite_difference(grid_n)
    exact_y = np.sin(np.pi * grid_x)
    error = np.max(np.abs(grid_y - exact_y))
    errors.append(error)

    if i == 0:
        print(f"{grid_n:8d} {error:20.6e} {'---':>12}")
    else:
        order = np.log(errors[i - 1] / error) / np.log(2)
        print(f"{grid_n:8d} {error:20.6e} {order:12.4f}")

# Plot numerical and exact solutions
plt.figure(figsize=(9, 5))
plt.plot(x, y_numerical, "o-", label="Numerical solution")
plt.plot(x, y_exact, "--", label="Exact solution")
plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Finite Difference Method")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("finite_difference_method.png", dpi=300)
plt.show()

# Plot absolute error
plt.figure(figsize=(9, 5))
plt.plot(x, absolute_error, "o-", label="Absolute error")
plt.xlabel("x")
plt.ylabel("Absolute error")
plt.title("Finite Difference Error")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("finite_difference_error.png", dpi=300)
plt.show()

# Plot convergence
plt.figure(figsize=(8, 5))
plt.loglog(grid_sizes, errors, "o-", label="Maximum error")
plt.loglog(
    grid_sizes,
    errors[0] * (np.array(grid_sizes[0]) / np.array(grid_sizes))**2,
    "--",
    label="Second-order reference"
)
plt.xlabel("Number of subintervals")
plt.ylabel("Maximum error")
plt.title("Finite Difference Convergence")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("finite_difference_convergence.png", dpi=300)
plt.show()
