import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# Number of interior grid points
grid_sizes = [20, 50, 100, 200, 400, 800]

# Exact first three eigenvalues
exact_values = np.array([
    np.pi**2,
    (2 * np.pi)**2,
    (3 * np.pi)**2
])

errors = np.zeros((len(grid_sizes), 3))
step_sizes = []

for i, N in enumerate(grid_sizes):

    h = 1 / (N + 1)
    step_sizes.append(h)

    # Finite difference matrix
    diagonal = np.full(N, 2 / h**2)
    off_diagonal = np.full(N - 1, -1 / h**2)

    # Calculate the first three numerical eigenvalues
    numerical_values, _ = eigh_tridiagonal(
        diagonal,
        off_diagonal,
        select="i",
        select_range=(0, 2)
    )

    # Calculate absolute errors
    errors[i, :] = np.abs(numerical_values - exact_values)

    print(f"\nGrid points: {N}")
    for k in range(3):
        print(
            f"Eigenvalue {k + 1}: "
            f"Numerical = {numerical_values[k]:.8f}, "
            f"Exact = {exact_values[k]:.8f}, "
            f"Error = {errors[i, k]:.8e}"
        )

# Plot errors
plt.figure(figsize=(9, 6))

for k in range(3):
    plt.loglog(
        step_sizes,
        errors[:, k],
        marker="o",
        label=f"Eigenvalue {k + 1}"
    )

plt.xlabel("Grid size h")
plt.ylabel("Absolute error")
plt.title("Accuracy of Numerical Eigenvalues")
plt.grid(True, which="both")
plt.legend()
plt.show()
