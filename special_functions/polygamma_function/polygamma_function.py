import numpy as np
import matplotlib.pyplot as plt
from scipy.special import polygamma, gamma
import math


# --------------------------------------------------
# Polygamma Function
#
# Definition:
# psi^(n)(x) = d^(n+1)/dx^(n+1) [ln Gamma(x)]
#
# n = 0, 1, 2, ...
# x > 0
# --------------------------------------------------


# --------------------------------------------------
# 1. Numerical Evaluation
# --------------------------------------------------

x_values = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]

orders = [0, 1, 2, 3, 4]

print("Numerical values of the Polygamma function")
print("-" * 60)

for n in orders:
    print(f"\nOrder n = {n}")
    for x in x_values:
        value = polygamma(n, x)
        print(f"psi^({n})({x}) = {value:.12f}")

# --------------------------------------------------
# 2. Verification of Special Values
# --------------------------------------------------

print("\n\nVerification of Special Values")
print("-" * 60)

gamma_constant = np.euler_gamma

special_values = [
    ("psi(1) = -gamma", polygamma(0, 1), -gamma_constant),
    (
        "psi(1/2) = -gamma - 2 ln(2)",
        polygamma(0, 0.5),
        -gamma_constant - 2 * math.log(2)
    ),
    ("psi^(1)(1) = pi^2/6", polygamma(1, 1), math.pi**2 / 6),
    (
        "psi^(1)(1/2) = pi^2/2",
        polygamma(1, 0.5),
        math.pi**2 / 2
    ),
]

for description, numerical, exact in special_values:
    error = abs(numerical - exact)
    print(f"\n{description}")
    print(f"Numerical value = {numerical:.12f}")
    print(f"Exact value     = {exact:.12f}")
    print(f"Absolute error  = {error:.3e}")

# --------------------------------------------------
# 3. Verification of Recurrence Relations
# --------------------------------------------------

print("\n\nVerification of Recurrence Relations")
print("-" * 60)

test_x_values = [0.5, 1.0, 1.5, 2.0, 3.0]
orders = [0, 1, 2, 3, 4]

for n in orders:
    print(f"\nOrder n = {n}")
    max_error = 0.0

    for x in test_x_values:
        lhs = polygamma(n, x + 1)

        rhs = (
            polygamma(n, x)
            + (-1)**n * math.factorial(n) / x**(n + 1)
        )

        error = abs(lhs - rhs)
        max_error = max(max_error, error)

        print(
            f"x = {x:3.1f} | "
            f"LHS = {lhs: .12f} | "
            f"RHS = {rhs: .12f} | "
            f"Error = {error:.3e}"
        )

    print(f"Maximum error for n = {n}: {max_error:.3e}")
# --------------------------------------------------
# 4. Series Representation and Convergence
# --------------------------------------------------

print("\n\nSeries Representation and Convergence")
print("-" * 60)

# Mathematical representation:
#
# psi^(n)(x) = (-1)^(n+1) n! *
#              sum_{k=0}^{infinity} 1/(x+k)^(n+1)
#
# Valid for n >= 1 and x > 0.


def polygamma_series(n, x, N):
    """Approximate the polygamma function using N terms."""
    partial_sum = sum(
        1.0 / (x + k)**(n + 1)
        for k in range(N)
    )

    return (-1)**(n + 1) * math.factorial(n) * partial_sum


x_test = 1.5
orders_to_test = [1, 2, 3]
N_values = [10, 50, 100, 500, 1000, 5000]

for n in orders_to_test:
    exact_value = polygamma(n, x_test)

    print(f"\nOrder n = {n}, x = {x_test}")
    print(f"Exact value = {exact_value:.12f}")
    print("-" * 60)
    print(f"{'N':>8} {'Approximation':>20} {'Absolute Error':>20}")

    for N in N_values:
        approximation = polygamma_series(n, x_test, N)
        error = abs(exact_value - approximation)

        print(
            f"{N:8d} "
            f"{approximation:20.12f} "
            f"{error:20.6e}"
        )
# --------------------------------------------------
# 5. Complete Monotonicity
# --------------------------------------------------

print("\n\nComplete Monotonicity Verification")
print("-" * 60)

# For n >= 1:
#
# (-1)^(n+1) psi^(n)(x) > 0, for x > 0.


x_grid = np.linspace(0.1, 10.0, 200)
orders_to_check = [1, 2, 3, 4, 5]

for n in orders_to_check:
    signed_values = (-1)**(n + 1) * polygamma(n, x_grid)

    minimum_value = np.min(signed_values)
    all_positive = np.all(signed_values > 0)

    print(f"\nOrder n = {n}")
    print(f"Minimum signed value = {minimum_value:.12e}")
    print(f"All values positive? {all_positive}")

print("\nComplete monotonicity sign pattern verified numerically.")


# --------------------------------------------------
# 6. Graphical Visualization
# --------------------------------------------------

print("\n\nGenerating Polygamma Function Plots...")

x_plot = np.linspace(0.1, 5.0, 500)
plot_orders = [0, 1, 2, 3]

fig, axes = plt.subplots(2, 2, figsize=(12, 9))

for n, ax in zip(plot_orders, axes.flat):
    y_plot = polygamma(n, x_plot)

    ax.plot(x_plot, y_plot, linewidth=2)
    ax.set_title(f"Polygamma Function: Order n = {n}")
    ax.set_xlabel("x")
    ax.set_ylabel(f"psi^({n})(x)")
    ax.grid(True, alpha=0.3)

fig.suptitle("Polygamma Functions of Different Orders", fontsize=15)
plt.tight_layout()

plt.savefig("polygamma_function/polygamma_orders.png", dpi=300)
plt.show()

print("Plot saved as polygamma_function/polygamma_orders.png")

# --------------------------------------------------
# 7. Logarithmic Comparison of Polygamma Functions
# --------------------------------------------------

print("\nGenerating logarithmic comparison plot...")

x_log = np.linspace(0.1, 5.0, 500)
plot_orders = [0, 1, 2, 3]

plt.figure(figsize=(10, 6))

for n in plot_orders:
    y = np.abs(polygamma(n, x_log))
    plt.semilogy(
        x_log,
        y,
        linewidth=2,
        label=f"Order n = {n}"
    )

plt.xlabel("x")
plt.ylabel(r"$|\psi^{(n)}(x)|$")
plt.title("Logarithmic Comparison of Polygamma Functions")
plt.grid(True, which="both", alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig(
    "polygamma_function/polygamma_log_comparison.png",
    dpi=300
)
plt.show()

print("Logarithmic comparison plot saved successfully.")

# --------------------------------------------------
# 8. Asymptotic Approximation and Error Analysis
# --------------------------------------------------

print("\n\nAsymptotic Approximation and Error Analysis")
print("-" * 70)

# Leading-order asymptotic formula:
#
# psi^(n)(x) ~ (-1)^(n-1) (n-1)! / x^n
#
# Valid for n >= 1 as x tends to infinity.


x_values_asymptotic = [5, 10, 20, 50, 100]
orders_asymptotic = [1, 2, 3, 4]

for n in orders_asymptotic:
    print(f"\nOrder n = {n}")
    print("-" * 70)

    print(
        f"{'x':>8} {'Exact value':>20} "
        f"{'Approximation':>20} {'Relative Error':>18}"
    )

    for x in x_values_asymptotic:
        exact = polygamma(n, x)

        approximation = (
            (-1)**(n - 1) * math.factorial(n - 1) / x**n
        )

        relative_error = abs(exact - approximation) / abs(exact)

        print(
            f"{x:8.1f} "
            f"{exact:20.12e} "
            f"{approximation:20.12e} "
            f"{relative_error:18.6e}"
        )
