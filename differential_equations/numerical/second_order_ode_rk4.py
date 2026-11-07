
import numpy as np
import matplotlib.pyplot as plt


# Homogeneous second-order ODE:
# y'' + 3y' + 2y = 0
#
# Initial conditions:
# y(0) = 1, y'(0) = 0
#
# Convert to a first-order system:
# y' = z
# z' = -3z - 2y


def system(x, Y):
    y, z = Y
    return np.array([z, -3*z - 2*y])


# Exact solutions
def exact_y(x):
    return 2*np.exp(-x) - np.exp(-2*x)


def exact_z(x):
    return -2*np.exp(-x) + 2*np.exp(-2*x)


# Classical fourth-order Runge-Kutta method
def rk4(f, x0, Y0, x_end, h):

    if h <= 0:
        raise ValueError("Step size must be positive.")

    if x_end <= x0:
        raise ValueError("The final point must exceed the initial point.")

    n = int(np.ceil((x_end - x0) / h))
    x = np.linspace(x0, x_end, n + 1)
    h_actual = (x_end - x0) / n

    Y = np.zeros((n + 1, len(Y0)))
    Y[0] = Y0

    for i in range(n):
        hi = x[i + 1] - x[i]

        k1 = f(x[i], Y[i])
        k2 = f(x[i] + hi/2, Y[i] + hi*k1/2)
        k3 = f(x[i] + hi/2, Y[i] + hi*k2/2)
        k4 = f(x[i] + hi, Y[i] + hi*k3)

        Y[i + 1] = Y[i] + (hi/6)*(k1 + 2*k2 + 2*k3 + k4)

    return x, Y, h_actual


# Initial-value problem
x0 = 0.0
x_end = 5.0
Y0 = np.array([1.0, 0.0])
h = 0.01

# Numerical solution
x, Y, h_actual = rk4(system, x0, Y0, x_end, h)

y_numerical = Y[:, 0]
z_numerical = Y[:, 1]

# Exact solutions
y_true = exact_y(x)
z_true = exact_z(x)

# Absolute errors
y_errors = np.abs(y_numerical - y_true)
z_errors = np.abs(z_numerical - z_true)


# Display numerical results
print("HOMOGENEOUS SECOND-ORDER ODE")
print("=" * 65)
print("Equation: y'' + 3y' + 2y = 0")
print("Initial conditions: y(0) = 1, y'(0) = 0")
print("Method: Classical fourth-order Runge-Kutta")
print(f"Requested step size: {h}")
print(f"Actual step size: {h_actual:.6f}")

print("\nNUMERICAL RESULTS")
print("-" * 65)
print(f"{'x':>8} {'Numerical y':>16} {'Exact y':>16} "
      f"{'Error in y':>14}")

for i in range(0, len(x), 50):
    print(
        f"{x[i]:8.2f} "
        f"{y_numerical[i]:16.10f} "
        f"{y_true[i]:16.10f} "
        f"{y_errors[i]:14.6e}"
    )

# Include the final point
i = len(x) - 1
print(
    f"{x[i]:8.2f} "
    f"{y_numerical[i]:16.10f} "
    f"{y_true[i]:16.10f} "
    f"{y_errors[i]:14.6e}"
)


# Error summary
print("\nERROR ANALYSIS")
print("-" * 40)
print(f"Maximum error in y:  {np.max(y_errors):.12e}")
print(f"Maximum error in y': {np.max(z_errors):.12e}")
print(f"Final error in y:    {y_errors[-1]:.12e}")
print(f"Final error in y':   {z_errors[-1]:.12e}")


# Convergence study
step_sizes = [0.2, 0.1, 0.05, 0.025]
max_errors = []

print("\nFOURTH-ORDER CONVERGENCE STUDY")
print("-" * 55)
print(f"{'Step size':>12} {'Maximum error':>20} {'Observed order':>18}")

for i, step in enumerate(step_sizes):

    x_test, Y_test, _ = rk4(system, x0, Y0, x_end, step)

    error = np.max(
        np.abs(Y_test[:, 0] - exact_y(x_test))
    )

    max_errors.append(error)

    if i == 0:
        print(f"{step:12.5f} {error:20.8e} {'---':>18}")
    else:
        order = np.log(max_errors[i-1] / error) / np.log(2)
        print(f"{step:12.5f} {error:20.8e} {order:18.4f}")


# Plot 1: Numerical and exact solutions
plt.figure(figsize=(9, 5))

plt.plot(x, y_numerical, label="Numerical solution (RK4)")
plt.plot(x, y_true, "--", label="Exact solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Homogeneous Second-Order ODE: RK4 Solution")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("second_order_ode_rk4.png", dpi=300)
plt.show()


# Plot 2: Absolute errors
plt.figure(figsize=(9, 5))

plt.semilogy(x, y_errors, label="Error in y")
plt.semilogy(x, z_errors, label="Error in y'")

plt.xlabel("x")
plt.ylabel("Absolute error (log scale)")
plt.title("RK4 Error Analysis")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("second_order_ode_rk4_error.png", dpi=300)
plt.show()


# Plot 3: Convergence study
plt.figure(figsize=(9, 5))

plt.loglog(step_sizes, max_errors, "o-", label="RK4 maximum error")

reference = [
    max_errors[0] * (step / step_sizes[0])**4
    for step in step_sizes
]

plt.loglog(step_sizes, reference, "--", label="Fourth-order reference")

plt.xlabel("Step size")
plt.ylabel("Maximum absolute error")
plt.title("Fourth-Order Convergence of RK4")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("second_order_ode_rk4_convergence.png", dpi=300)
plt.show()
