
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# Define the interval and discretization
L = 1.0
N = 100
h = L / N

# Construct the finite difference matrix for -y''
diagonal = 2 * np.ones(N - 1) / h**2
off_diagonal = -np.ones(N - 2) / h**2

# Compute the first three eigenvalues and eigenvectors
eigenvalues, eigenvectors = eigh_tridiagonal(
    diagonal,
    off_diagonal,
    select="i",
    select_range=(0, 2)
)

# Numerical grid
x = np.linspace(0, L, N + 1)

print("Numerical Eigenvalue Problem")
print("=" * 65)

# Compare numerical and exact eigenvalues
for n, eigenvalue in enumerate(eigenvalues, start=1):

    exact = (n * np.pi)**2
    absolute_error = abs(eigenvalue - exact)
    relative_error = absolute_error / exact

    print(f"\nEigenvalue {n}")
    print(f"Numerical value: {eigenvalue:.10f}")
    print(f"Exact value:     {exact:.10f}")
    print(f"Absolute error:  {absolute_error:.8e}")
    print(f"Relative error:  {relative_error:.8e}")

# Plot the numerical eigenfunctions
plt.figure(figsize=(9, 6))

for i in range(3):

    eigenfunction = np.zeros(N + 1)
    eigenfunction[1:N] = eigenvectors[:, i]

    # Normalize the eigenfunction
    eigenfunction /= np.max(np.abs(eigenfunction))

    # Choose a consistent positive orientation
    max_index = np.argmax(np.abs(eigenfunction))

    if eigenfunction[max_index] < 0:
        eigenfunction *= -1

    plt.plot(
        x,
        eigenfunction,
        linewidth=2,
        label=f"Eigenfunction {i + 1}"
    )

plt.xlabel("x")
plt.ylabel("Eigenfunction")
plt.title("Numerical Eigenfunctions")
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save the figure
plt.savefig("eigenvalue_problem.png", dpi=300)

plt.show()
