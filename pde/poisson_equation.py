
import numpy as np
import matplotlib.pyplot as plt
import os


# --------------------------------------------------
# Two-Dimensional Poisson Equation
#
# u_xx + u_yy = f(x,y)
#
# Exact solution:
# u(x,y) = sin(pi*x) * sin(pi*y)
#
# Source term:
# f(x,y) = -2*pi^2*sin(pi*x)*sin(pi*y)
#
# Boundary conditions:
# u = 0 on all four boundaries
# --------------------------------------------------


# Physical domain
Lx = 1.0
Ly = 1.0

# Numerical parameters
Nx = 30
Ny = 30

tolerance = 1e-8
max_iterations = 50000


def exact_solution(X, Y):
    """Return the exact analytical solution."""
    return np.sin(np.pi * X) * np.sin(np.pi * Y)


def source_term(X, Y):
    """Return the source term of the Poisson equation."""
    return -2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)


def gauss_seidel(Nx, Ny, tolerance=1e-8, max_iterations=50000):

    x = np.linspace(0.0, Lx, Nx + 1)
    y = np.linspace(0.0, Ly, Ny + 1)

    dx = Lx / Nx
    dy = Ly / Ny

    X, Y = np.meshgrid(x, y)

    f = source_term(X, Y)

    # Initial approximation
    u = np.zeros((Ny + 1, Nx + 1))

    # Gauss-Seidel iteration
    for iteration in range(1, max_iterations + 1):

        max_change = 0.0

        for j in range(1, Ny):

            for i in range(1, Nx):

                old_value = u[j, i]

                u[j, i] = (
                    dy**2 * (u[j, i + 1] + u[j, i - 1])
                    + dx**2 * (u[j + 1, i] + u[j - 1, i])
                    - dx**2 * dy**2 * f[j, i]
                ) / (2.0 * (dx**2 + dy**2))

                change = abs(u[j, i] - old_value)

                max_change = max(max_change, change)

        if max_change < tolerance:
            break

    return x, y, X, Y, u, iteration, max_change


# --------------------------------------------------
# Numerical Solution
# --------------------------------------------------

x, y, X, Y, u, iterations, final_change = gauss_seidel(
    Nx, Ny, tolerance, max_iterations
)

u_exact = exact_solution(X, Y)

error = np.abs(u - u_exact)

maximum_error = np.max(error)
rmse = np.sqrt(np.mean(error**2))


# --------------------------------------------------
# Numerical Results
# --------------------------------------------------

print("Gauss-Seidel Numerical Results")
print("=" * 50)

print(f"Grid size: {Nx} x {Ny}")
print(f"Number of iterations: {iterations}")
print(f"Final maximum change: {final_change:.10e}")
print(f"Maximum absolute error: {maximum_error:.10e}")
print(f"Root mean square error: {rmse:.10e}")


# --------------------------------------------------
# Visualization of Numerical and Exact Solutions
# --------------------------------------------------

os.makedirs("images", exist_ok=True)

# 1. Numerical solution: 3D surface
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, u, cmap="viridis", edgecolor="none")

ax.set_title("Numerical Solution of the 2D Poisson Equation")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u(x,y)")

plt.tight_layout()
plt.savefig("images/poisson_equation_2d_numerical.png", dpi=300)
plt.close()


# 2. Exact solution: 3D surface
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, u_exact, cmap="viridis", edgecolor="none")

ax.set_title("Exact Solution of the 2D Poisson Equation")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u(x,y)")

plt.tight_layout()
plt.savefig("images/poisson_equation_2d_exact.png", dpi=300)
plt.close()


# 3. Contour plot of numerical solution
plt.figure(figsize=(8, 6))

contour = plt.contourf(X, Y, u, levels=30, cmap="viridis")
plt.colorbar(contour, label="u(x,y)")

plt.title("Numerical Solution: Contour Plot")
plt.xlabel("x")
plt.ylabel("y")

plt.tight_layout()
plt.savefig("images/poisson_equation_2d_contour.png", dpi=300)
plt.close()


# 4. Absolute error distribution
plt.figure(figsize=(8, 6))

contour = plt.contourf(X, Y, error, levels=30, cmap="inferno")
plt.colorbar(contour, label="Absolute Error")

plt.title("Absolute Error Distribution")
plt.xlabel("x")
plt.ylabel("y")

plt.tight_layout()
plt.savefig("images/poisson_equation_2d_error.png", dpi=300)
plt.close()


# --------------------------------------------------
# Grid Convergence Study
# --------------------------------------------------

grid_sizes = [10, 20, 30, 40]

maximum_errors = []
rms_errors = []

print("\nGrid Convergence Study")
print("=" * 65)
print(f"{'Grid':<12}{'Maximum Error':<25}{'RMSE':<20}")

for N in grid_sizes:

    _, _, X_grid, Y_grid, u_grid, _, _ = gauss_seidel(
        N, N, tolerance, max_iterations
    )

    exact_grid = exact_solution(X_grid, Y_grid)

    error_grid = np.abs(u_grid - exact_grid)

    max_error = np.max(error_grid)
    grid_rmse = np.sqrt(np.mean(error_grid**2))

    maximum_errors.append(max_error)
    rms_errors.append(grid_rmse)

    print(
        f"{N} x {N:<8}"
        f"{max_error:<25.10e}"
        f"{grid_rmse:<20.10e}"
    )


# Calculate observed convergence orders
orders = []

for k in range(1, len(grid_sizes)):

    order = np.log(
        maximum_errors[k - 1] / maximum_errors[k]
    ) / np.log(
        grid_sizes[k] / grid_sizes[k - 1]
    )

    orders.append(order)

print("\nObserved Convergence Orders:")

for k, order in enumerate(orders):

    print(
        f"{grid_sizes[k]} -> {grid_sizes[k + 1]}: "
        f"{order:.4f}"
    )


# 5. Log-log convergence plot
dx_values = [1.0 / N for N in grid_sizes]

plt.figure(figsize=(8, 6))

plt.loglog(
    dx_values,
    maximum_errors,
    "o-",
    label="Maximum absolute error"
)

plt.loglog(
    dx_values,
    rms_errors,
    "s--",
    label="RMSE"
)

plt.xlabel("Grid spacing")
plt.ylabel("Numerical error")
plt.title("Grid Convergence Study: 2D Poisson Equation")

plt.grid(True, which="both", linestyle="--")
plt.legend()

plt.tight_layout()
plt.savefig("images/poisson_equation_2d_convergence.png", dpi=300)
plt.close()

print("\nAll five figures saved successfully in the images folder.")
