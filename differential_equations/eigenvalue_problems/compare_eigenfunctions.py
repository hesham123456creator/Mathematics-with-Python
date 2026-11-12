
"""
Compare numerical and exact eigenfunctions for -y'' = lambda*y
with homogeneous Dirichlet boundary conditions on [0, 1].
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# Number of interior grid points
N = 100

# Grid
h = 1 / (N + 1)
x = np.linspace(h, 1 - h, N)
x_full = np.linspace(0, 1, N + 2)

# Construct the finite difference matrix
diagonal = np.full(N, 2 / h**2)
off_diagonal = np.full(N - 1, -1 / h**2)

# Compute the first three eigenvalues and eigenvectors
eigenvalues, eigenvectors = eigh_tridiagonal(
    diagonal,
    off_diagonal,
    select="i",
    select_range=(0, 2)
)

print("Comparison of Numerical and Exact Eigenfunctions")
print("=" * 60)

plt.figure(figsize=(10, 6))

for k in range(3):

    n = k + 1

    # Numerical eigenfunction
    numerical = np.zeros(N + 2)
    numerical[1:-1] = eigenvectors[:, k]

    # Normalize and choose positive orientation
    numerical /= np.max(np.abs(numerical))

    if numerical[1] < 0:
        numerical *= -1

       # Exact eigenfunction
    exact = np.sin(n * np.pi * x_full)
    exact /= np.max(np.abs(exact))

    # Calculate maximum absolute error
    error = np.max(np.abs(numerical - exact))

    print(f"\nEigenfunction {n}")
    print(f"Numerical eigenvalue: {eigenvalues[k]:.8f}")
    print(f"Exact eigenvalue:     {(n * np.pi)**2:.8f}")
    print(f"Maximum absolute error: {error:.8e}")

    # Combined comparison plot
    plt.plot(
        x_full,
        numerical,
        "--",
        label=f"Numerical eigenfunction {n}"
    )

    plt.plot(
        x_full,
        exact,
        label=f"Exact eigenfunction {n}"
    )

plt.xlabel("x")
plt.ylabel("Eigenfunction")
plt.title("Numerical vs Exact Eigenfunctions")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig(
    "eigenfunction_comparison.png",
    dpi=300
)

plt.show()
plt.close()
