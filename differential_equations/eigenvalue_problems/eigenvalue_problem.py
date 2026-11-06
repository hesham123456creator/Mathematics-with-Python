import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# Define the interval
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

print("Numerical and exact eigenvalues:")

for n, eigenvalue in enumerate(eigenvalues, start=1):
    exact = (n * np.pi)**2
    print(
        f"n = {n}, Numerical = {eigenvalue:.8f}, "
        f"Exact = {exact:.8f}"
    )

# Plot the eigenfunctions
for i in range(3):
    eigenfunction = np.zeros(N + 1)
    eigenfunction[1:N] = eigenvectors[:, i]

    # Normalize the eigenfunction
    eigenfunction /= np.max(np.abs(eigenfunction))

    # Make the first eigenfunction positive near the center
    if eigenfunction[N // 2] < 0:
        eigenfunction *= -1

    plt.plot(x, eigenfunction, label=f"Eigenfunction {i + 1}")

plt.xlabel("x")
plt.ylabel("Eigenfunction")
plt.title("Eigenfunctions of a Boundary Value Problem")
plt.legend()
plt.grid()
plt.show()
