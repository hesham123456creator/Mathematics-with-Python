
import os
import time
import numpy as np
import matplotlib.pyplot as plt


# --------------------------------------------------
# Comparison: Jacobi vs Gauss-Seidel
# 2D Poisson Equation
# --------------------------------------------------

Lx = 1.0
Ly = 1.0

Nx = 30
Ny = 30

tolerance = 1e-8
max_iterations = 50000


def exact_solution(X, Y):
    return np.sin(np.pi * X) * np.sin(np.pi * Y)


def source_term(X, Y):
    return -2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)


def jacobi(Nx, Ny, tolerance, max_iterations):

    x = np.linspace(0.0, Lx, Nx + 1)
    y = np.linspace(0.0, Ly, Ny + 1)

    dx = Lx / Nx
    dy = Ly / Ny

    X, Y = np.meshgrid(x, y)
    f = source_term(X, Y)

    u = np.zeros((Ny + 1, Nx + 1))
    history = []

    for iteration in range(1, max_iterations + 1):

        u_new = u.copy()

        u_new[1:Ny, 1:Nx] = (
            dy**2 * (u[1:Ny, 2:Nx+1] + u[1:Ny, 0:Nx-1])
            + dx**2 * (u[2:Ny+1, 1:Nx] + u[0:Ny-1, 1:Nx])
            - dx**2 * dy**2 * f[1:Ny, 1:Nx]
        ) / (2.0 * (dx**2 + dy**2))

        max_change = np.max(np.abs(u_new - u))
        history.append(max_change)

        u = u_new

        if max_change < tolerance:
            break

    return x, y, X, Y, u, iteration, max_change, history


def gauss_seidel(Nx, Ny, tolerance, max_iterations):

    x = np.linspace(0.0, Lx, Nx + 1)
    y = np.linspace(0.0, Ly, Ny + 1)

    dx = Lx / Nx
    dy = Ly / Ny

    X, Y = np.meshgrid(x, y)
    f = source_term(X, Y)

    u = np.zeros((Ny + 1, Nx + 1))
    history = []

    for iteration in range(1, max_iterations + 1):

        max_change = 0.0

        for j in range(1, Ny):
            for i in range(1, Nx):

                old_value = u[j, i]

                u[j, i] = (
                    dy**2 * (u[j, i+1] + u[j, i-1])
                    + dx**2 * (u[j+1, i] + u[j-1, i])
                    - dx**2 * dy**2 * f[j, i]
                ) / (2.0 * (dx**2 + dy**2))

                change = abs(u[j, i] - old_value)
                max_change = max(max_change, change)

        history.append(max_change)

        if max_change < tolerance:
            break

    return x, y, X, Y, u, iteration, max_change, history


def calculate_errors(u, X, Y):

    u_exact = exact_solution(X, Y)
    error = np.abs(u - u_exact)

    maximum_error = np.max(error)
    rmse = np.sqrt(np.mean(error**2))

    return maximum_error, rmse


def main():

    os.makedirs("images", exist_ok=True)

    print("2D Poisson Equation: Solver Comparison")
    print("=" * 65)

    print(f"Grid size: {Nx} x {Ny}")
    print(f"Tolerance: {tolerance:.1e}")

    # Warm-up both methods before timing
    jacobi(Nx, Ny, tolerance, max_iterations)
    gauss_seidel(Nx, Ny, tolerance, max_iterations)

    results = {}

    # Run each method three times and use median execution time
    for name, solver in [
        ("Jacobi", jacobi),
        ("Gauss-Seidel", gauss_seidel)
    ]:

        times = []

        for _ in range(3):

            start = time.perf_counter()

            result = solver(
                Nx, Ny, tolerance, max_iterations
            )

            elapsed = time.perf_counter() - start
            times.append(elapsed)

        x, y, X, Y, u, iterations, final_change, history = result

        maximum_error, rmse = calculate_errors(u, X, Y)

        results[name] = {
            "iterations": iterations,
            "time": np.median(times),
            "maximum_error": maximum_error,
            "rmse": rmse,
            "final_change": final_change,
            "history": history,
            "u": u
        }

    # --------------------------------------------------
    # Numerical comparison table
    # --------------------------------------------------

    print("\nNumerical Results")
    print("=" * 85)

    print(
        f"{'Method':<16}"
        f"{'Iterations':<14}"
        f"{'Time (s)':<14}"
        f"{'Maximum Error':<22}"
        f"{'RMSE':<18}"
    )

    for name, result in results.items():

        print(
            f"{name:<16}"
            f"{result['iterations']:<14}"
            f"{result['time']:<14.6f}"
            f"{result['maximum_error']:<22.10e}"
            f"{result['rmse']:<18.10e}"
        )

    # --------------------------------------------------
    # Performance comparison
    # --------------------------------------------------

    jacobi_iterations = results["Jacobi"]["iterations"]
    gs_iterations = results["Gauss-Seidel"]["iterations"]

    reduction = (
        (jacobi_iterations - gs_iterations)
        / jacobi_iterations * 100
    )

    speedup = (
        results["Jacobi"]["time"]
        / results["Gauss-Seidel"]["time"]
    )

    print("\nPerformance Analysis")
    print("=" * 50)

    print(f"Iteration reduction: {reduction:.2f}%")
    print(f"Measured time ratio (Jacobi / GS): {speedup:.3f}")

    # --------------------------------------------------
    # 1. Iteration comparison
    # --------------------------------------------------

    methods = list(results.keys())
    iterations = [results[m]["iterations"] for m in methods]

    plt.figure(figsize=(8, 6))
    plt.bar(methods, iterations)

    plt.ylabel("Number of Iterations")
    plt.title("Jacobi vs Gauss-Seidel: Iterations")

    plt.tight_layout()
    plt.savefig("images/poisson_iterations_comparison.png", dpi=300)
    plt.close()

    # --------------------------------------------------
    # 2. Execution time comparison
    # --------------------------------------------------

    times = [results[m]["time"] for m in methods]

    plt.figure(figsize=(8, 6))
    plt.bar(methods, times)

    plt.ylabel("Execution Time (seconds)")
    plt.title("Jacobi vs Gauss-Seidel: Execution Time")

    plt.tight_layout()
    plt.savefig("images/poisson_time_comparison.png", dpi=300)
    plt.close()

    # --------------------------------------------------
    # 3. Numerical accuracy comparison
    # --------------------------------------------------

    max_errors = [results[m]["maximum_error"] for m in methods]
    rms_errors = [results[m]["rmse"] for m in methods]

    positions = np.arange(len(methods))
    width = 0.35

    plt.figure(figsize=(9, 6))

    plt.bar(
        positions - width/2,
        max_errors,
        width,
        label="Maximum Error"
    )

    plt.bar(
        positions + width/2,
        rms_errors,
        width,
        label="RMSE"
    )

    plt.xticks(positions, methods)
    plt.ylabel("Error")
    plt.title("Numerical Accuracy Comparison")
    plt.legend()

    plt.tight_layout()
    plt.savefig("images/poisson_accuracy_comparison.png", dpi=300)
    plt.close()

    # --------------------------------------------------
    # 4. Convergence history
    # --------------------------------------------------

    plt.figure(figsize=(9, 6))

    for name, result in results.items():

        plt.semilogy(
            range(1, len(result["history"]) + 1),
            result["history"],
            label=name
        )

    plt.axhline(
        tolerance,
        linestyle="--",
        label="Tolerance"
    )

    plt.xlabel("Iteration")
    plt.ylabel("Maximum Change")
    plt.title("Convergence History: Jacobi vs Gauss-Seidel")

    plt.grid(True, which="both", linestyle="--")
    plt.legend()

    plt.tight_layout()
    plt.savefig("images/poisson_convergence_comparison.png", dpi=300)
    plt.close()

    print("\nFour comparison figures saved successfully.")
    print("Location: images/")


if __name__ == "__main__":
    main()
