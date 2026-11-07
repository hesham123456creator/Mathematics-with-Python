
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# Exact first three eigenvalues
exact_values = np.array([
    np.pi**2,
    (2 * np.pi)**2,
    (3 * np.pi)**2
])

# Different numbers of subintervals
grid_sizes = [20, 40, 80, 160, 320]

errors = np.zeros((len(grid_sizes), 3))

print("Numerical Eigenvalue Accuracy Analysis")
print("=" * 65)

for j, N in enumerate(grid_sizes):

    # Grid spacing
    h = 1.0 / N

    # Finite difference matrix for -y''
    diagonal = 2 * np.ones(N - 1) / h**2
    off_diagonal = -np.ones(N - 2) / h**2

    # Compute the first three eigenvalues
    numerical_values, _ = eigh_tridiagonal(
        diagonal,
        off_diagonal,
        select="i",
        select_range=(0, 2)
    )

    # Calculate absolute errors
    errors[j, :] = np.abs(numerical_values - exact_values)

    print(f"\nNumber of subintervals: {N}")

    for i in range(3):
        print(
            f"Eigenvalue {i + 1}: "
            f"Numerical = {numerical_values[i]:.10f}, "
            f"Exact = {exact_values[i]:.10f}, "
            f"Error = {errors[j, i]:.8e}"
        )

# Calculate observed convergence orders
print("\nObserved Convergence Orders")
print("=" * 65)

for i in range(3):
    print(f"\nEigenvalue {i + 1}:")

    for j in range(1, len(grid_sizes)):
        order = np.log(
            errors[j - 1, i] / errors[j, i]
        ) / np.log(
            grid_sizes[j] / grid_sizes[j - 1]
        )

        print(
            f"N = {grid_sizes[j - 1]} to {grid_sizes[j]}: "
            f"Order = {order:.4f}"
        )

# Plot absolute errors
plt.figure(figsize=(9, 6))

for i in range(3):
    plt.loglog(
        grid_sizes,
        errors[:, i],
        "o-",
        label=f"Eigenvalue {i + 1}"
    )

plt.xlabel("Number of subintervals (N)")
plt.ylabel("Absolute error")
plt.title("Convergence of Numerical Eigenvalues")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("eigenvalue_accuracy.png", dpi=300)
plt.show()

