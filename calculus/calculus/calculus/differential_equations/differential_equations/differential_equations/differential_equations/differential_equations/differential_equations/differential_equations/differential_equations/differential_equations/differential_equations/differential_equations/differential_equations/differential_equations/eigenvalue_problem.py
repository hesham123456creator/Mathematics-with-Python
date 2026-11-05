
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# Number of interior grid points
N = 200

# Step size
h = 1 / (N + 1)

# Construct the finite difference matrix
diagonal = np.full(N, 2 / h**2)
off_diagonal = np.full(N - 1, -1 / h**2)

# Calculate eigenvalues and eigenvectors
eigenvalues, eigenvectors = eigh_tridiagonal(
    diagonal,
    off_diagonal
)

# First five eigenvalues
print("Numerical and Exact Eigenvalues")
print("--------------------------------")

for n in range(5):

    numerical = eigenvalues[n]
    exact = ((n + 1) * np.pi)**2

    error = abs(numerical - exact)

    print(
        f"n={n+1}: "
        f"Numerical={numerical:.8f}, "
        f"Exact={exact:.8f}, "
        f"Error={error:.8e}"
    )

# Grid including boundary points
x = np.linspace(0, 1, N + 2)

# Plot first three eigenfunctions
plt.figure(figsize=(9, 5))

for n in range(3):

    y = np.zeros(N + 2)

    # Normalize eigenfunction
    y[1:-1] = eigenvectors[:, n]

    # Normalize to maximum absolute value 1
    y = y / np.max(np.abs(y))

    # Choose positive orientation
    if y[1] < 0:
        y = -y

    plt.plot(
        x,
        y,
        label=f"Eigenfunction {n+1}"
    )

plt.xlabel("x")
plt.ylabel("Eigenfunction")
plt.title("First Three Numerical Eigenfunctions")
plt.legend()
plt.grid(True)
plt.show()
