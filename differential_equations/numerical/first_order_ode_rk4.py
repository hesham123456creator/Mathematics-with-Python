
import numpy as np
import matplotlib.pyplot as plt


# Differential equation: y' = -2y
def f(x, y):
    return -2 * y


# Exact solution
def exact_solution(x):
    return np.exp(-2 * x)


# Classical fourth-order Runge-Kutta method
def rk4(f, x0, y0, x_end, h):
    if h <= 0:
        raise ValueError("Step size must be positive.")

    if x_end <= x0:
        raise ValueError("x_end must be greater than x0.")

    n = int(np.ceil((x_end - x0) / h))
    x = np.linspace(x0, x_end, n + 1)
    h_actual = x[1] - x[0]

    y = np.zeros(n + 1)
    y[0] = y0

    for i in range(n):
        k1 = f(x[i], y[i])
        k2 = f(x[i] + h_actual / 2, y[i] + h_actual * k1 / 2)
        k3 = f(x[i] + h_actual / 2, y[i] + h_actual * k2 / 2)
        k4 = f(x[i] + h_actual, y[i] + h_actual * k3)

        y[i + 1] = y[i] + (h_actual / 6) * (
            k1 + 2 * k2 + 2 * k3 + k4
        )

    return x, y


# Initial-value problem
x0 = 0.0
y0 = 1.0
x_end = 5.0
h = 0.01


# Compute numerical solution
x, y_numerical = rk4(f, x0, y0, x_end, h)

# Exact solution
y_exact = exact_solution(x)

# Absolute errors
errors = np.abs(y_numerical - y_exact)

print("FIRST-ORDER ODE: CLASSICAL RK4 METHOD")
print("=" * 55)
print("Differential equation: y' = -2y")
print("Initial condition: y(0) = 1")
print(f"Interval: [{x0}, {x_end}]")
print(f"Step size requested: {h}")
print(f"Number of steps: {len(x) - 1}")
print()

# Display selected numerical results
print(f"{'x':>8} {'Numerical RK4':>18} "
      f"{'Exact':>18} {'Absolute Error':>18}")
print("-" * 65)

indices = np.unique(
    np.append(np.arange(0, len(x), 50), len(x) - 1)
)

for i in indices:
    print(
        f"{x[i]:8.2f} "
        f"{y_numerical[i]:18.10f} "
        f"{y_exact[i]:18.10f} "
        f"{errors[i]:18.6e}"
    )

print("\nERROR ANALYSIS")
print("-" * 40)
print(f"Maximum absolute error: {np.max(errors):.12e}")
print(f"Final absolute error:   {errors[-1]:.12e}")
print(f"Final numerical value:  {y_numerical[-1]:.12e}")
print(f"Final exact value:      {y_exact[-1]:.12e}")


# Convergence study
step_sizes = [0.2, 0.1, 0.05, 0.025]
maximum_errors = []

for step in step_sizes:
    x_test, y_test = rk4(f, x0, y0, x_end, step)
    y_exact_test = exact_solution(x_test)

    maximum_errors.append(
        np.max(np.abs(y_test - y_exact_test))
    )

print("\nCONVERGENCE STUDY")
print("-" * 55)
print(f"{'Step size':>12} {'Maximum error':>20} {'Observed order':>18}")

for i, step in enumerate(step_sizes):
    if i == 0:
        order_text = "---"
    else:
        order = np.log(
            maximum_errors[i - 1] / maximum_errors[i]
        ) / np.log(step_sizes[i - 1] / step)
        order_text = f"{order:.4f}"

    print(
        f"{step:12.4f} "
        f"{maximum_errors[i]:20.10e} "
        f"{order_text:>18}"
    )


# Plot 1: Numerical and exact solutions
plt.figure(figsize=(9, 5))

plt.plot(x, y_exact, "--", label="Exact solution", linewidth=2)
plt.plot(x, y_numerical, label="Numerical RK4 solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("First-Order ODE: Classical RK4 Method")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("first_order_ode_rk4.png", dpi=300)
plt.show()


# Plot 2: Absolute error
plt.figure(figsize=(9, 5))

plt.semilogy(x, errors, label="RK4 absolute error")

plt.xlabel("x")
plt.ylabel("Absolute error (log scale)")
plt.title("RK4 Numerical Error")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("first_order_ode_rk4_error.png", dpi=300)
plt.show()


# Plot 3: Convergence behavior
plt.figure(figsize=(9, 5))

plt.loglog(
    step_sizes,
    maximum_errors,
    "o-",
    label="RK4 maximum error"
)

plt.xlabel("Step size")
plt.ylabel("Maximum absolute error")
plt.title("Convergence Study of the RK4 Method")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("first_order_ode_rk4_convergence.png", dpi=300)
plt.show()

