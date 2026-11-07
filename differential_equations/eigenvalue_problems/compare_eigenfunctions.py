

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# Problem interval and discretization
L = 1.0
N = 100
h = L / N

# Finite difference matrix for -y''
diagonal = 2 * np.ones(N - 1) / h**2
off_diagonal = -np.ones(N - 2) / h**2

# Compute the first three eigenvalues and eigenvectors
eigenvalues, eigenvectors = eigh_tridiagonal(
    diagonal,
    off_diagonal,
    select="i",
    select_range=(0, 2)
)

# Grid including boundary points
x = np.linspace(0, L, N + 1)

print("Comparison of Numerical and Exact Eigenfunctions")
print("=" * 60)

for i in range(3):

    n = i + 1

    # Numerical eigenfunction
    numerical = np.zeros(N + 1)
    numerical[1:N] = eigenvectors[:, i]

    # Normalize numerical eigenfunction
    numerical /= np.max(np.abs(numerical))

    # Choose positive orientation
    max_index = np.argmax(np.abs(numerical))

    if numerical[max_index] < 0:
        numerical *= -1

    # Exact eigenfunction
    exact = np.sin(n * np.pi * x)

    # Calculate maximum absolute error
    error = np.max(np.abs(numerical - exact))

    print(f"\nEigenfunction {n}")
    print(f"Maximum absolute error: {error:.8e}")

    # Plot comparison
    plt.figure(figsize=(8, 5))

    plt.plot(
        x,
        numerical,
        label="Numerical eigenfunction",
        linewidth=2
    )

    plt.plot(
        x,
        exact,
        "--",
        label="Exact eigenfunction",
        linewidth=2
    )

    plt.xlabel("x")
    plt.ylabel("Eigenfunction")
    plt.title(f"Comparison of Eigenfunction {n}")

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Save the figure
    plt.savefig(
        f"eigenfunction_comparison_{n}.png",
        dpi=300
    )

    plt.show()
    plt.close()
