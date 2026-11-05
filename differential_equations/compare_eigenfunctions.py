import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# Problem interval
L = 1.0
N = 100
h = L / N

# Finite difference matrix for -y''
diagonal = 2 * np.ones(N - 1) / h**2
off_diagonal = -np.ones(N - 2) / h**2

# Compute the first three eigenpairs
eigenvalues, eigenvectors = eigh_tridiagonal(
    diagonal,
    off_diagonal,
    select="i",
    select_range=(0, 2)
)

# Grid including boundary points
x = np.linspace(0, L, N + 1)

# Compare numerical and exact eigenfunctions
for i in range(3):

    n = i + 1

    # Numerical eigenfunction
    numerical = np.zeros(N + 1)
    numerical[1:N] = eigenvectors[:, i]

    # Normalize numerical eigenfunction
    numerical /= np.max(np.abs(numerical))

    # Correct the sign
    if numerical[N // 2] < 0:
        numerical *= -1

    # Exact eigenfunction
    exact = np.sin(n * np.pi * x)

    # Calculate maximum error
    error = np.max(np.abs(numerical - exact))

    print(f"Eigenfunction {n}: Maximum error = {error:.8e}")

    # Plot comparison
    plt.figure()
    plt.plot(x, numerical, label="Numerical eigenfunction")
    plt.plot(x, exact, "--", label="Exact eigenfunction")

    plt.xlabel("x")
    plt.ylabel("Eigenfunction")
    plt.title(f"Comparison of Eigenfunction {n}")
    plt.legend()
    plt.grid()
    plt.show()
