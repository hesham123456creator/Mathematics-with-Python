
import numpy as np
import matplotlib.pyplot as plt
import time

# Grid parameters
Nx = 30
Ny = 30

Lx = 1.0
Ly = 1.0

dx = Lx / Nx
dy = Ly / Ny

x = np.linspace(0, Lx, Nx + 1)

# Boundary conditions
def initialize_grid():
    u = np.zeros((Ny + 1, Nx + 1))
    u[-1, :] = np.sin(np.pi * x)
    return u

# Jacobi method
def jacobi_method(tolerance=1e-6, max_iterations=20000):

    u = initialize_grid()

    for iteration in range(max_iterations):

        u_new = u.copy()

        u_new[1:-1, 1:-1] = (
            dy**2 * (
                u[1:-1, 2:] + u[1:-1, :-2]
            )
            + dx**2 * (
                u[2:, 1:-1] + u[:-2, 1:-1]
            )
        ) / (2 * (dx**2 + dy**2))

        error = np.max(np.abs(u_new - u))
        u = u_new

        if error < tolerance:
            break

    return u, iteration + 1

# Gauss-Seidel method
def gauss_seidel_method(
    tolerance=1e-6,
    max_iterations=20000
):

    u = initialize_grid()

    for iteration in range(max_iterations):

        max_change = 0.0

        for j in range(1, Ny):
            for i in range(1, Nx):

                old_value = u[j, i]

                u[j, i] = (
                    dy**2 * (
                        u[j, i+1] + u[j, i-1]
                    )
                    + dx**2 * (
                        u[j+1, i] + u[j-1, i]
                    )
                ) / (2 * (dx**2 + dy**2))

                change = abs(u[j, i] - old_value)
                max_change = max(max_change, change)

        if max_change < tolerance:
            break

    return u, iteration + 1

# Run Jacobi method
start = time.perf_counter()
u_jacobi, iter_jacobi = jacobi_method()
time_jacobi = time.perf_counter() - start

# Run Gauss-Seidel method
start = time.perf_counter()
u_gs, iter_gs = gauss_seidel_method()
time_gs = time.perf_counter() - start

# Compare results
difference = np.max(np.abs(u_jacobi - u_gs))

print("Jacobi Method")
print("Iterations:", iter_jacobi)
print("Execution time:", time_jacobi)

print("\nGauss-Seidel Method")
print("Iterations:", iter_gs)
print("Execution time:", time_gs)

print("\nMaximum difference:", difference)

# Plot both numerical solutions
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].imshow(
    u_jacobi,
    origin="lower",
    extent=[0, Lx, 0, Ly],
    aspect="auto"
)
axes[0].set_title("Jacobi Method")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")

axes[1].imshow(
    u_gs,
    origin="lower",
    extent=[0, Lx, 0, Ly],
    aspect="auto"
)
axes[1].set_title("Gauss-Seidel Method")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")

plt.tight_layout()
plt.show()
