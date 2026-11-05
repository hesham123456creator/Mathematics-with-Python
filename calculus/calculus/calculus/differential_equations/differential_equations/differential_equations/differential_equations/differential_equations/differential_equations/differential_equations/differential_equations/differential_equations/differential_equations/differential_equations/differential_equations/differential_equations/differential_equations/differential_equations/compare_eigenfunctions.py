import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# Number of interior grid points
N = 100

# Grid
h = 1 / (N + 1)
x = np.linspace(h, 1 - h, N)

# Construct the finite difference matrix
diagonal = np.full(N, 2 / h**2)
off_diagonal = np.full(N - 1, -1 / h**2)

# Calculate the first three eigenvalues and eigenvectors
eigenvalues, eigenvectors = eigh_tridiagonal(
    diagonal,
    off_diagonal,
    select="i",
    select_range=(0, 2)
)

# Add boundary points
x_full = np.linspace(0, 1, N + 2)

plt.figure(figsize=(10, 6))

for k in range(3):

    # Numerical eigenfunction
    numerical = np.zeros(N + 2)
    numerical[1:-1] = eigenvectors[:, k]

    # Normalize and choose positive orientation
    numerical = numerical / np.max(np.abs(numerical))
    if numerical[1] < 0:
        numerical = -numerical

    # Exact eigenfunction
    exact = np.sin((k + 1) * np.pi * x_full)

    plt.plot(
        x_full,
        numerical,
        "--",
        label=f"Numerical eigenfunction {k + 1}"
    )

    plt.plot(
        x_full,
        exact,
        label=f"Exact eigenfunction {k + 1}"
    )

    print(f"Eigenvalue {k + 1}: {eigenvalues[k]:.8f}")

plt.xlabel("x")
plt.ylabel("Eigenfunction")
plt.title("Numerical vs Exact Eigenfunctions")
plt.grid(True)
plt.legend()
plt.show()
