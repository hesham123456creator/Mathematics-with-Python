
import numpy as np
import matplotlib.pyplot as plt

# ==================================================
# Two-Dimensional Laplace Equation
# Gauss-Seidel Finite Difference Method
# ==================================================

# Domain parameters
Lx = 1.0
Ly = 1.0

Nx = 30
Ny = 30

dx = Lx / Nx
dy = Ly / Ny

# Create grid
x = np.linspace(0, Lx, Nx + 1)
y = np.linspace(0, Ly, Ny + 1)
X, Y = np.meshgrid(x, y)

# Exact analytical solution
def exact_solution(X, Y):
    return (
        np.sin(np.pi * X)
        * np.sinh(np.pi * Y)
        / np.sinh(np.pi)
    )


# Gauss-Seidel numerical solver
def gauss_seidel(Nx, Ny, tolerance=1e-8,
                 max_iterations=50000):

    dx = Lx / Nx
    dy = Ly / Ny

    x = np.linspace(0, Lx, Nx + 1)
    y = np.linspace(0, Ly, Ny + 1)

    u = np.zeros((Ny + 1, Nx + 1))

    # Boundary conditions
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    u[0, :] = 0.0
    u[-1, :] = np.sin(np.pi * x)

    # Iterative process
    for iteration in range(max_iterations):

        max_change = 0.0

        for j in range(1, Ny):
            for i in range(1, Nx):

                old_value = u[j, i]

                u[j, i] = (
                    dy**2 * (u[j, i+1] + u[j, i-1])
                    + dx**2 * (u[j+1, i] + u[j-1, i])
                ) / (2 * (dx**2 + dy**2))

                change = abs(u[j, i] - old_value)
                max_change = max(max_change, change)

        if max_change < tolerance:
            break

    return x, y, u, iteration + 1, max_change


# Solve the problem
x, y, u, iterations, final_change = gauss_seidel(Nx, Ny)

X, Y = np.meshgrid(x, y)

# Analytical solution
u_exact = exact_solution(X, Y)

# Error analysis
error = np.abs(u - u_exact)

maximum_error = np.max(error)
rmse = np.sqrt(np.mean(error**2))

# Numerical results
print("\nGauss-Seidel Numerical Results")
print("=" * 45)

print("Grid size:", Nx, "x", Ny)
print("Number of iterations:", iterations)
print(f"Final maximum change: {final_change:.10e}")
print(f"Maximum absolute error: {maximum_error:.10e}")
print(f"Root mean square error: {rmse:.10e}")


# ==================================================
# Visualization of Numerical and Exact Solutions
# ==================================================

# Create output directory
import os
os.makedirs("images", exist_ok=True)

# 1. Numerical solution: 3D surface
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, u, cmap="viridis")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("Temperature")
ax.set_title("Numerical Solution: Laplace Equation")

plt.tight_layout()
plt.savefig("images/laplace_equation_2d_numerical.png",
            dpi=300)
plt.show()


# 2. Analytical solution: 3D surface
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, u_exact, cmap="viridis")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("Temperature")
ax.set_title("Exact Analytical Solution")

plt.tight_layout()
plt.savefig("images/laplace_equation_2d_exact.png",
            dpi=300)
plt.show()


# 3. Numerical solution: contour plot
plt.figure(figsize=(8, 6))

plt.contourf(X, Y, u, levels=25, cmap="viridis")
plt.colorbar(label="Temperature")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Numerical Temperature Distribution")

plt.tight_layout()
plt.savefig("images/laplace_equation_2d_contour.png",
            dpi=300)
plt.show()


# 4. Absolute error distribution
plt.figure(figsize=(8, 6))

plt.contourf(X, Y, error, levels=25, cmap="inferno")
plt.colorbar(label="Absolute Error")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Absolute Error: Numerical vs Exact")

plt.tight_layout()
plt.savefig("images/laplace_equation_2d_error.png",
            dpi=300)
plt.show()


# ==================================================
# Grid Convergence Study
# ==================================================

grid_sizes = [10, 20, 30, 40]

maximum_errors = []
rms_errors = []

print("\nGrid Convergence Study")
print("=" * 65)
print(f"{'Grid':<12}{'Maximum Error':<22}{'RMSE':<22}")

for N in grid_sizes:

    x_grid, y_grid, u_grid, iters, change = (
        gauss_seidel(N, N)
    )

    X_grid, Y_grid = np.meshgrid(x_grid, y_grid)

    exact_grid = exact_solution(X_grid, Y_grid)

    error_grid = np.abs(u_grid - exact_grid)

    max_error = np.max(error_grid)
    rms_error = np.sqrt(np.mean(error_grid**2))

    maximum_errors.append(max_error)
    rms_errors.append(rms_error)

    print(
        f"{N} x {N:<8}"
        f"{max_error:<22.10e}"
        f"{rms_error:<22.10e}"
    )


# Estimate convergence order
orders = []

for k in range(len(grid_sizes) - 1):

    order = np.log(
        maximum_errors[k] / maximum_errors[k + 1]
    ) / np.log(
        grid_sizes[k + 1] / grid_sizes[k]
    )

    orders.append(order)

print("\nObserved convergence orders:")

for k, order in enumerate(orders):
    print(
        f"{grid_sizes[k]} -> {grid_sizes[k+1]}: "
        f"{order:.4f}"
    )


# 5. Convergence plot
plt.figure(figsize=(8, 6))

plt.loglog(
    grid_sizes,
    maximum_errors,
    "o-",
    label="Maximum absolute error"
)

plt.loglog(
    grid_sizes,
    rms_errors,
    "s--",
    label="RMSE"
)

plt.xlabel("Number of grid intervals")
plt.ylabel("Error")
plt.title("Grid Convergence Study")
plt.legend()
plt.grid(True, which="both", linestyle="--")

plt.tight_layout()
plt.savefig("images/laplace_equation_2d_convergence.png",
            dpi=300)
plt.show()

print("\nAll figures saved successfully in the images folder.")
