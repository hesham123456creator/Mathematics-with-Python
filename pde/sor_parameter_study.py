
import numpy as np
import matplotlib.pyplot as plt
import time

from sor_poisson import sor_poisson


def main():
    omega_values = np.arange(1.70, 1.91, 0.01)

    iterations_list = []
    times_list = []
    errors_list = []

    print("SOR Relaxation Parameter Study")
    print("=" * 75)
    print(
        f"{'Omega':<12}"
        f"{'Iterations':<15}"
        f"{'Time (s)':<15}"
        f"{'Maximum Error':<20}"
    )

    for omega in omega_values:
        results = sor_poisson(
            n=30,
            omega=omega,
            tol=1e-8,
            max_iter=50000
        )

        iterations = results[5]
        max_error = results[7]
        elapsed_time = results[9]

        iterations_list.append(iterations)
        times_list.append(elapsed_time)
        errors_list.append(max_error)

        print(
            f"{omega:<12.2f}"
            f"{iterations:<15}"
            f"{elapsed_time:<15.6f}"
            f"{max_error:<20.10e}"
        )

    best_index = np.argmin(iterations_list)
    best_omega = omega_values[best_index]

    print("\nOptimal Parameter (by iteration count)")
    print("=" * 50)
    print(f"Best omega: {best_omega:.2f}")
    print(f"Minimum iterations: {iterations_list[best_index]}")
    print(f"Maximum error: {errors_list[best_index]:.10e}")

    plt.figure(figsize=(8, 6))
    plt.plot(
        omega_values,
        iterations_list,
        marker="o",
        linewidth=2
    )
    plt.xlabel("Relaxation parameter (omega)")
    plt.ylabel("Number of iterations")
    plt.title("SOR Convergence versus Relaxation Parameter")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        "images/sor_refined_iterations.png",
        dpi=300
    )
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.plot(
        omega_values,
        times_list,
        marker="s",
        linewidth=2
    )
    plt.xlabel("Relaxation parameter (omega)")
    plt.ylabel("Execution time (seconds)")
    plt.title("SOR Execution Time versus Relaxation Parameter")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        "images/sor_refined_time.png",
        dpi=300
    )
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.semilogy(
        omega_values,
        errors_list,
        marker="^",
        linewidth=2
    )
    plt.xlabel("Relaxation parameter (omega)")
    plt.ylabel("Maximum absolute error")
    plt.title("SOR Accuracy versus Relaxation Parameter")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        "images/sor_refined_accuracy.png",
        dpi=300
    )
    plt.close()

    print("\nThree parameter-study figures saved successfully.")


if __name__ == "__main__":
    main()
