import numpy as np
import matplotlib.pyplot as plt
import time
import os


# --------------------------------------------------
# Two-Dimensional Poisson Equation
#
# -Delta u = f(x,y)
#
# Exact solution:
# u(x,y) = sin(pi*x) * sin(pi*y)
#
# Source term:
# f(x,y) = 2*pi^2*sin(pi*x)*sin(pi*y)
#
# Boundary conditions:
# u = 0 on all four boundaries
# --------------------------------------------------


# Physical domain
Lx = 1.0
Ly = 1.0

# Number of intervals
Nx = 30
Ny = 30

# Numerical parameters
tolerance = 1e-8
max_iterations = 50000

# SOR relaxation parameter
omega = 1.5

# --------------------------------------------------
# Exact Solution and Source Term
# --------------------------------------------------

def exact_solution(X, Y):
    """Return the exact analytical solution."""
    return np.sin(np.pi * X) * np.sin(np.pi * Y)


def source_term(X, Y):
    """Return the source term."""
    return 2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)


# --------------------------------------------------
# Common Grid and Initial Approximation
# --------------------------------------------------

def initialize_problem(Nx, Ny):

    x = np.linspace(0.0, Lx, Nx + 1)
    y = np.linspace(0.0, Ly, Ny + 1)

    dx = Lx / Nx
    dy = Ly / Ny

    X, Y = np.meshgrid(x, y)

    F = source_term(X, Y)

    # Initial approximation
    U = np.zeros((Ny + 1, Nx + 1))

    return x, y, X, Y, dx, dy, F, U

# --------------------------------------------------
# Relative Residual Norm
# --------------------------------------------------

def relative_residual(U, F, dx, dy):

    residual = (
        (U[1:-1, 2:] - 2.0 * U[1:-1, 1:-1] + U[1:-1, :-2]) / dx**2
        +
        (U[2:, 1:-1] - 2.0 * U[1:-1, 1:-1] + U[:-2, 1:-1]) / dy**2
        + F[1:-1, 1:-1]
    )

    numerator = np.linalg.norm(residual.ravel())
    denominator = np.linalg.norm(F[1:-1, 1:-1].ravel())

    return numerator / denominator

# --------------------------------------------------
# Jacobi Method
# --------------------------------------------------

def jacobi_method():

    x, y, X, Y, dx, dy, F, U = initialize_problem(Nx, Ny)

    start_time = time.perf_counter()

    residual = relative_residual(U, F, dx, dy)
    residual_history = [residual]

    for iteration in range(1, max_iterations + 1):

        U_new = U.copy()

        U_new[1:-1, 1:-1] = (
            dy**2 * (U[1:-1, 2:] + U[1:-1, :-2])
            + dx**2 * (U[2:, 1:-1] + U[:-2, 1:-1])
            + dx**2 * dy**2 * F[1:-1, 1:-1]
        ) / (2.0 * (dx**2 + dy**2))

        U = U_new

        residual = relative_residual(U, F, dx, dy)
        residual_history.append(residual)

        if residual < tolerance:
            break

    elapsed_time = time.perf_counter() - start_time

    return U, iteration, residual, elapsed_time, residual_history

# --------------------------------------------------
# Gauss-Seidel Method
# --------------------------------------------------

def gauss_seidel_method():

    x, y, X, Y, dx, dy, F, U = initialize_problem(Nx, Ny)

    start_time = time.perf_counter()

    residual = relative_residual(U, F, dx, dy)
    residual_history = [residual]

    for iteration in range(1, max_iterations + 1):

        for j in range(1, Ny):
            for i in range(1, Nx):

                U[j, i] = (
                    dy**2 * (U[j, i + 1] + U[j, i - 1])
                    + dx**2 * (U[j + 1, i] + U[j - 1, i])
                    + dx**2 * dy**2 * F[j, i]
                ) / (2.0 * (dx**2 + dy**2))

        residual = relative_residual(U, F, dx, dy)
        residual_history.append(residual)

        if residual < tolerance:
            break

    elapsed_time = time.perf_counter() - start_time

    return U, iteration, residual, elapsed_time, residual_history

# --------------------------------------------------
# Successive Over-Relaxation (SOR) Method
# --------------------------------------------------

def sor_method():

    x, y, X, Y, dx, dy, F, U = initialize_problem(Nx, Ny)

    start_time = time.perf_counter()

    residual = relative_residual(U, F, dx, dy)
    residual_history = [residual]

    for iteration in range(1, max_iterations + 1):

        for j in range(1, Ny):
            for i in range(1, Nx):

                old_value = U[j, i]

                gs_value = (
                    dy**2 * (U[j, i + 1] + U[j, i - 1])
                    + dx**2 * (U[j + 1, i] + U[j - 1, i])
                    + dx**2 * dy**2 * F[j, i]
                ) / (2.0 * (dx**2 + dy**2))

                U[j, i] = (
                    (1.0 - omega) * old_value
                    + omega * gs_value
                )

        residual = relative_residual(U, F, dx, dy)
        residual_history.append(residual)

        if residual < tolerance:
            break

    elapsed_time = time.perf_counter() - start_time

    return U, iteration, residual, elapsed_time, residual_history

# --------------------------------------------------
# Conjugate Gradient Method
# --------------------------------------------------

def conjugate_gradient_method():

    x, y, X, Y, dx, dy, F, U = initialize_problem(Nx, Ny)

    start_time = time.perf_counter()

    # Right-hand side: interior points only
    b = F[1:-1, 1:-1].copy()

    # Initial residual
    r = b.copy()

    # Initial search direction
    p = r.copy()

    # Squared residual norm
    rr_old = np.sum(r * r)

    b_norm = np.linalg.norm(b)

    residual = np.sqrt(rr_old) / b_norm
    residual_history = [residual]

    for iteration in range(1, max_iterations + 1):

        # Apply the discrete positive-definite operator A
        # A approximates -Laplacian

        p_padded = np.pad(p, 1, mode="constant")

        Ap = (
            (2.0 * p - p_padded[1:-1, 2:] - p_padded[1:-1, :-2]) / dx**2
            +
            (2.0 * p - p_padded[2:, 1:-1] - p_padded[:-2, 1:-1]) / dy**2
        )

        # Step length
        alpha = rr_old / np.sum(p * Ap)

        # Update solution
        U[1:-1, 1:-1] += alpha * p

        # Update residual
        r -= alpha * Ap

        rr_new = np.sum(r * r)

        # Recalculate the true relative residual
        residual = relative_residual(U, F, dx, dy)
        residual_history.append(residual)

        if residual < tolerance:
            break


        # Update search direction
        beta = rr_new / rr_old

        p = r + beta * p

        rr_old = rr_new

    elapsed_time = time.perf_counter() - start_time

    return U, iteration, residual, elapsed_time, residual_history

# --------------------------------------------------
# Main Comparison
# --------------------------------------------------

def main():

    x, y, X, Y, dx, dy, F, U_initial = initialize_problem(Nx, Ny)

    u_exact = exact_solution(X, Y)

    methods = {
        "Jacobi": jacobi_method,
        "Gauss-Seidel": gauss_seidel_method,
        "SOR": sor_method,
        "Conjugate Gradient": conjugate_gradient_method
    }

    results = {}

    print("\nComparison of Four Numerical Methods")
    print("=" * 95)
    print(
        f"{'Method':<24}"
        f"{'Iterations':<15}"
        f"{'Residual':<18}"
        f"{'Maximum Error':<20}"
        f"{'Time (s)':<12}"
    )
    print("-" * 95)

    for name, method in methods.items():

        U, iterations, residual, elapsed_time, history = method()

        error = np.abs(U - u_exact)
        max_error = np.max(error)
        rmse = np.sqrt(np.mean(error**2))

        results[name] = {
            "solution": U,
            "iterations": iterations,
            "residual": residual,
            "time": elapsed_time,
            "max_error": max_error,
            "rmse": rmse,
            "history": history
        }

        print(
            f"{name:<24}"
            f"{iterations:<15}"
            f"{residual:<18.6e}"
            f"{max_error:<20.6e}"
            f"{elapsed_time:<12.6f}"
        )

    print("\nComparison completed successfully.")

    # --------------------------------------------------
    # Convergence History
    # --------------------------------------------------

    plt.figure(figsize=(10, 6))

    for name, data in results.items():
        plt.semilogy(
            data["history"],
            label=name
        )

    plt.xlabel("Iteration")
    plt.ylabel("Relative Residual")
    plt.title("Convergence Comparison of Four Methods")

    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    os.makedirs("images", exist_ok=True)

    plt.savefig(
        "images/poisson_methods_convergence.png",
        dpi=300
    )

    plt.show()
    
    # --------------------------------------------------
    # Numerical Solution Comparison
    # --------------------------------------------------

    fig = plt.figure(figsize=(14, 10))

    for index, (name, data) in enumerate(results.items(), start=1):

        ax = fig.add_subplot(2, 2, index, projection="3d")

        surface = ax.plot_surface(
            X,
            Y,
            data["solution"],
            cmap="viridis",
            edgecolor="none"
        )

        ax.set_title(name)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("u(x,y)")

    fig.suptitle("Numerical Solutions of the 2D Poisson Equation")

    plt.tight_layout()

    plt.savefig(
        "images/poisson_methods_solutions.png",
        dpi=300
    )

    plt.show()

    # --------------------------------------------------
    # Absolute Error Comparison
    # --------------------------------------------------

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for ax, (name, data) in zip(axes.flat, results.items()):

        error = np.abs(data["solution"] - u_exact)

        contour = ax.contourf(
            X,
            Y,
            error,
            levels=20,
            cmap="inferno"
        )

        fig.colorbar(contour, ax=ax, label="Absolute Error")

        ax.set_title(name)
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    fig.suptitle("Absolute Error Distribution")

    plt.tight_layout()

    plt.savefig(
        "images/poisson_methods_errors.png",
        dpi=300
    )

    plt.show()

if __name__ == "__main__":
    main()
