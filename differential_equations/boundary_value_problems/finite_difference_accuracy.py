
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

# Exact solution
def exact_solution(x):
    return np.sin(np.pi * x)


# Finite Difference Method
def finite_difference(N):

    h = 1.0 / N
    x = np.linspace(0.0, 1.0, N + 1)

    # Construct the tridiagonal matrix
    diagonal = -2.0 * np.ones(N - 1)
    off_diagonal = np.ones(N - 2)

    ab = np.zeros((3, N - 1))
    ab[0, 1:] = off_diagonal
    ab[1, :] = diagonal
    ab[2, :-1] = off_diagonal

    # Right-hand side
    rhs = -h**2 * np.pi**2 * np.sin(np.pi * x[1:N])

    # Solve the linear system
    y_interior = solve_banded((1, 1), ab, rhs)

    # Include boundary conditions
    y = np.zeros(N + 1)
    y[1:N] = y_interior

    return x, y


# --------------------------------------------------
# 1. Test different grid sizes
# --------------------------------------------------

grid_sizes = [10, 20, 40, 80, 160, 320]

errors = []

print("Finite Difference Accuracy Analysis")
print("=" * 65)

for N in grid_sizes:

    x, y_numerical = finite_difference(N)
    y_exact = exact_solution(x)

    error = np.max(np.abs(y_numerical - y_exact))
    errors.append(error)

    print(
        f"N = {N:3d}, "
        f"Maximum error = {error:.10e}"
    )


# --------------------------------------------------
# 2. Calculate observed convergence orders
# --------------------------------------------------

print("\nObserved Convergence Orders")
print("=" * 65)

orders = []

for i in range(1, len(grid_sizes)):

    order = np.log(errors[i - 1] / errors[i]) / np.log(
        grid_sizes[i] / grid_sizes[i - 1]
    )

    orders.append(order)

    print(
        f"N = {grid_sizes[i - 1]:3d} to {grid_sizes[i]:3d}: "
        f"Order = {order:.6f}"
    )


# --------------------------------------------------
# 3. Plot error and second-order reference
# --------------------------------------------------

reference = [
    errors[0] * (grid_sizes[0] / N)**2
    for N in grid_sizes
]

plt.figure(figsize=(9, 6))

plt.loglog(
    grid_sizes,
    errors,
    "o-",
    label="Finite Difference Error",
    linewidth=2
)

plt.loglog(
    grid_sizes,
    reference,
    "--",
    label="Second-order reference",
    linewidth=2
)

plt.xlabel("Number of subintervals (N)")
plt.ylabel("Maximum absolute error")
plt.title("Convergence of the Finite Difference Method")

plt.legend()
plt.grid(True, which="both")
plt.tight_layout()

plt.savefig("finite_difference_accuracy.png", dpi=300)

plt.show()
plt.close()

