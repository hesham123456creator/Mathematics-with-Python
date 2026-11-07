
import numpy as np
import matplotlib.pyplot as plt


# Nonhomogeneous second-order differential equation:
# y'' + 3y' + 2y = exp(-x)
# Initial conditions: y(0) = 0, y'(0) = 0


# Convert the second-order equation into a first-order system:
# y' = z
# z' = exp(-x) - 3z - 2y

def system(x, Y):
    y, z = Y

    dydx = z
    dzdx = np.exp(-x) - 3 * z - 2 * y

    return np.array([dydx, dzdx])


# Exact analytical solution
def exact_solution(x):
    return (x - 1) * np.exp(-x) + np.exp(-2 * x)


def exact_derivative(x):
    return (2 - x) * np.exp(-x) - 2 * np.exp(-2 * x)


# Classical fourth-order Runge-Kutta method
def rk4(f, x0, Y0, x_end, h):
    if h <= 0:
        raise ValueError("Step size must be positive.")

    if x_end <= x0:
        raise ValueError("x_end must be greater than x0.")

    n = int(np.ceil((x_end - x0) / h))
    x = np.linspace(x0, x_end, n + 1)
    h_actual = x[1] - x[0]

    Y = np.zeros((n + 1, len(Y0)))
    Y[0] = Y0

    for i in range(n):
        k1 = f(x[i], Y[i])
        k2 = f(x[i] + h_actual / 2, Y[i] + h_actual * k1 / 2)
        k3 = f(x[i] + h_actual / 2, Y[i] + h_actual * k2 / 2)
        k4 = f(x[i] + h_actual, Y[i] + h_actual * k3)

        Y[i + 1] = Y[i] + (h_actual / 6) * (
            k1 + 2 * k2 + 2 * k3 + k4
        )

    return x, Y


# Initial-value problem
x0 = 0.0
Y0 = np.array([0.0, 0.0])

x_end = 5.0
h = 0.01


# Numerical solution
x, Y = rk4(system, x0, Y0, x_end, h)

y_numerical = Y[:, 0]
z_numerical = Y[:, 1]


# Exact solutions
y_exact = exact_solution(x)
z_exact = exact_derivative(x)


# Absolute errors
error_y = np.abs(y_numerical - y_exact)
error_z = np.abs(z_numerical - z_exact)


# Display selected numerical results
print("NONHOMOGENEOUS SECOND-ORDER ODE")
print("=" * 75)
print("Equation: y'' + 3y' + 2y = exp(-x)")
print("Initial conditions: y(0) = 0, y'(0) = 0")
print("Method: Classical fourth-order Runge-Kutta")
print(f"Step size: {h}")
print(f"Number of steps: {len(x) - 1}")

print("\nNUMERICAL RESULTS")
print("-" * 75)
print(f"{'x':>7} {'Numerical y':>14} {'Exact y':>14} "
      f"{'Numerical y prime':>19} {'Exact y prime':>15}")

indices = np.unique(
    np.append(np.arange(0, len(x), 50), len(x) - 1)
)

for i in indices:
    print(
        f"{x[i]:7.2f} "
        f"{y_numerical[i]:14.8f} "
        f"{y_exact[i]:14.8f} "
        f"{z_numerical[i]:19.8f} "
        f"{z_exact[i]:15.8f}"
    )


# Error analysis
print("\nERROR ANALYSIS")
print("-" * 45)
print(f"Maximum error in y:  {np.max(error_y):.12e}")
print(f"Maximum error in y': {np.max(error_z):.12e}")
print(f"Final error in y:    {error_y[-1]:.12e}")
print(f"Final error in y':   {error_z[-1]:.12e}")


# Convergence study
step_sizes = [0.2, 0.1, 0.05, 0.025]
maximum_errors = []

for step in step_sizes:
    x_test, Y_test = rk4(system, x0, Y0, x_end, step)
    y_test = Y_test[:, 0]
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
plt.title("Nonhomogeneous Second-Order ODE: RK4 Method")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("nonhomogeneous_ode_rk4.png", dpi=300)
plt.show()


# Plot 2: Absolute errors
plt.figure(figsize=(9, 5))

plt.semilogy(x, error_y, label="Error in y")
plt.semilogy(x, error_z, label="Error in y'")

plt.xlabel("x")
plt.ylabel("Absolute error (log scale)")
plt.title("RK4 Error Analysis")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("nonhomogeneous_ode_rk4_error.png", dpi=300)
plt.show()


# Plot 3: Convergence study
plt.figure(figsize=(9, 5))

plt.loglog(
    step_sizes,
    maximum_errors,
    "o-",
    label="RK4 maximum error"
)

plt.xlabel("Step size")
plt.ylabel("Maximum absolute error")
plt.title("RK4 Convergence for a Nonhomogeneous ODE")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("nonhomogeneous_ode_rk4_convergence.png", dpi=300)
plt.show()
