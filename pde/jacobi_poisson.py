
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
# Boundary conditions:
# u = 0 on all four boundaries
# --------------------------------------------------

Lx = 1.0
Ly = 1.0

Nx = 30
Ny = 30

tolerance = 1e-8
max_iterations = 50000


def exact_solution(X, Y):
    """Return the exact analytical solution."""
    return np.sin(np.pi * X) * np.sin(np.pi * Y)


def source_term(X, Y):
    """Return the source term."""
    return -2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)


def jacobi(Nx, Ny, tolerance=1e-8, max_iterations=50000):
    """Solve the 2D Poisson equation using Jacobi iteration."""

    x = np.linspace(0.0, Lx, Nx + 1)
    y = np.linspace(0.0, Ly, Ny + 1)

    dx = Lx / Nx
    dy = Ly / Ny

    X, Y = np.meshgrid(x, y)
    f = source_term(X, Y)

    # Initial approximation
    u = np.zeros((Ny + 1, Nx + 1))

    # Jacobi iteration
    for iteration in range(1, max_iterations + 1):

        u_new = u.copy()

        u_new[1:Ny, 1:Nx] = (
            dy**2 * (u[1:Ny, 2:Nx+1] + u[1:Ny, 0:Nx-1])
            + dx**2 * (u[2:Ny+1, 1:Nx] + u[0:Ny-1, 1:Nx])
            - dx**2 * dy**2 * f[1:Ny, 1:Nx]
        ) / (2.0 * (dx**2 + dy**2))

        max_change = np.max(np.abs(u_new - u))

        u = u_new

        if max_change < tolerance:
            break

    return x, y, X, Y, u, iteration, max_change


def main():

    x, y, X, Y, u, iterations, final_change = jacobi(
        Nx, Ny, tolerance, max_iterations
    )

    u_exact = exact_solution(X, Y)
    error = np.abs(u - u_exact)

    maximum_error = np.max(error)
    rmse = np.sqrt(np.mean(error**2))

    print("Jacobi Numerical Results")
    print("=" * 50)

    print(f"Grid size: {Nx} x {Ny}")
    print(f"Number of iterations: {iterations}")
    print(f"Final maximum change: {final_change:.10e}")
    print(f"Maximum absolute error: {maximum_error:.10e}")
    print(f"Root mean square error: {rmse:.10e}")

    os.makedirs("images", exist_ok=True)

    # Numerical solution
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_surface(X, Y, u, cmap="viridis", edgecolor="none")

    ax.set_title("Jacobi Solution of the 2D Poisson Equation")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("u(x,y)")

    plt.tight_layout()
    plt.savefig("images/poisson_jacobi_2d.png", dpi=300)
    plt.close()

    # Absolute error
    plt.figure(figsize=(8, 6))

    contour = plt.contourf(X, Y, error, levels=30, cmap="inferno")
    plt.colorbar(contour, label="Absolute Error")

    plt.title("Jacobi Absolute Error Distribution")
    plt.xlabel("x")
    plt.ylabel("y")

    plt.tight_layout()
    plt.savefig("images/poisson_jacobi_error.png", dpi=300)
    plt.close()

    print("\nTwo figures saved successfully in the images folder.")


if __name__ == "__main__":
    main()
