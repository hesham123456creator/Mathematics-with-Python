import numpy as np
import matplotlib.pyplot as plt

# Exact solution
def exact_solution(x):
    return np.sin(np.pi * x)

# Solve the BVP using the finite difference method
def finite_difference(N):
    h = 1 / N
    x = np.linspace(0, 1, N + 1)

    A = np.zeros((N - 1, N - 1))
    rhs = np.zeros(N - 1)

    for i in range(N - 1):
        A[i, i] = -2

        if i > 0:
            A[i, i - 1] = 1

        if i < N - 2:
            A[i, i + 1] = 1

        rhs[i] = -h**2 * np.pi**2 * np.sin(np.pi * x[i + 1])

    y = np.zeros(N + 1)
    y[1:N] = np.linalg.solve(A, rhs)

    return x, y

# Test different grid sizes
grid_sizes = [10, 20, 40, 80, 160]
errors = []

for N in grid_sizes:
    x, y_numerical = finite_difference(N)
    y_exact = exact_solution(x)

    error = np.max(np.abs(y_numerical - y_exact))
    errors.append(error)

    print(f"N = {N:3d}, Maximum error = {error:.8e}")

# Plot error against grid size
plt.loglog(grid_sizes, errors, "o-", label="Maximum error")
plt.loglog(
    grid_sizes,
    [errors[0] * (grid_sizes[0] / N)**2 for N in grid_sizes],
    "--",
    label="Second-order reference"
)

plt.xlabel("Number of subintervals (N)")
plt.ylabel("Maximum absolute error")
plt.title("Accuracy of the Finite Difference Method")
plt.legend()
plt.grid(True, which="both")
plt.show()
