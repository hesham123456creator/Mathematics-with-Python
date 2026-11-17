
import numpy as np
import matplotlib.pyplot as plt
import time


def exact_solution(x, y):
    return np.sin(np.pi * x) * np.sin(np.pi * y)


def sor_poisson(n=30, omega=1.5, tol=1e-8, max_iter=50000):
    """Solve the 2D Poisson equation using SOR."""

    h = 1.0 / (n + 1)
    x = np.linspace(0, 1, n + 2)
    y = np.linspace(0, 1, n + 2)
    X, Y = np.meshgrid(x, y)

    U = np.zeros((n + 2, n + 2))

    F = 2 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)

    start_time = time.perf_counter()

    for iteration in range(1, max_iter + 1):
        max_change = 0.0

        for i in range(1, n + 1):
            for j in range(1, n + 1):

                old_value = U[i, j]

                gs_value = 0.25 * (
                    U[i + 1, j]
                    + U[i - 1, j]
                    + U[i, j + 1]
                    + U[i, j - 1]
                    + h**2 * F[i, j]
                )

                U[i, j] = (
                    (1 - omega) * old_value
                    + omega * gs_value
                )

                change = abs(U[i, j] - old_value)
                max_change = max(max_change, change)

        if max_change < tol:
            break

    elapsed_time = time.perf_counter() - start_time

    U_exact = exact_solution(X, Y)
    error = np.abs(U - U_exact)

    max_error = np.max(error)
    rmse = np.sqrt(np.mean(error**2))

    return (
        X, Y, U, U_exact, error,
        iteration, max_change, max_error,
        rmse, elapsed_time
    )


def main():
    n = 30
    omega = 1.5
    tol = 1e-8

    results = sor_poisson(n=n, omega=omega, tol=tol)

    (
        X, Y, U, U_exact, error,
        iterations, final_change, max_error,
        rmse, elapsed_time
    ) = results

    print("SOR Numerical Results")
    print("=" * 50)
    print(f"Grid size: {n} x {n}")
    print(f"Relaxation parameter: {omega}")
    print(f"Number of iterations: {iterations}")
    print(f"Final maximum change: {final_change:.10e}")
    print(f"Maximum absolute error: {max_error:.10e}")
    print(f"Root mean square error: {rmse:.10e}")
    print(f"Execution time: {elapsed_time:.6f} seconds")

    plt.figure(figsize=(8, 6))
    contour = plt.contourf(X, Y, U, levels=30, cmap="viridis")
    plt.colorbar(contour, label="Numerical solution")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"SOR Solution (omega={omega})")
    plt.tight_layout()
    plt.savefig("images/poisson_sor_2d.png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 6))
    contour = plt.contourf(X, Y, error, levels=30, cmap="inferno")
    plt.colorbar(contour, label="Absolute error")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("SOR Absolute Error")
    plt.tight_layout()
    plt.savefig("images/poisson_sor_error.png", dpi=300)
    plt.close()

    print("\nTwo figures saved successfully.")


if __name__ == "__main__":
    main()
