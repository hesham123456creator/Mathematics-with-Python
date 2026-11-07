
import numpy as np
import matplotlib.pyplot as plt

# Differential equation: y' = y - x^2 + 1
def f(x, y):
    return y - x**2 + 1


# Exact analytical solution
def exact_solution(x):
    return (x + 1)**2 - 0.5 * np.exp(x)


# Initial-value problem
x0 = 0.0
y0 = 0.5
x_end = 2.0
n_steps = 10

# Discretize the interval
x = np.linspace(x0, x_end, n_steps + 1)
h = x[1] - x[0]


# Euler's method
y_euler = np.zeros(n_steps + 1)
y_euler[0] = y0

for i in range(n_steps):
    y_euler[i + 1] = (
        y_euler[i] + h * f(x[i], y_euler[i])
    )


# Fourth-order Runge-Kutta method
y_rk4 = np.zeros(n_steps + 1)
y_rk4[0] = y0

for i in range(n_steps):
    k1 = f(x[i], y_rk4[i])
    k2 = f(x[i] + h / 2, y_rk4[i] + h * k1 / 2)
    k3 = f(x[i] + h / 2, y_rk4[i] + h * k2 / 2)
    k4 = f(x[i] + h, y_rk4[i] + h * k3)

    y_rk4[i + 1] = y_rk4[i] + (h / 6) * (
        k1 + 2 * k2 + 2 * k3 + k4
    )


# Exact solution values
y_exact = exact_solution(x)


# Absolute errors
error_euler = np.abs(y_euler - y_exact)
error_rk4 = np.abs(y_rk4 - y_exact)


# Display numerical results
print("Numerical Solution of y' = y - x^2 + 1")
print("=" * 85)
print(f"{'x':>8} {'Euler':>16} {'RK4':>16} "
      f"{'Exact':>16} {'Euler Error':>14} {'RK4 Error':>14}")
print("-" * 85)

for i in range(n_steps + 1):
    print(
        f"{x[i]:8.2f} "
        f"{y_euler[i]:16.8f} "
        f"{y_rk4[i]:16.8f} "
        f"{y_exact[i]:16.8f} "
        f"{error_euler[i]:14.6e} "
        f"{error_rk4[i]:14.6e}"
    )


# Error summary
print("\nError Analysis")
print("-" * 40)
print(f"Step size: {h:.4f}")
print(f"Maximum Euler error: {np.max(error_euler):.10e}")
print(f"Maximum RK4 error:   {np.max(error_rk4):.10e}")
print(f"Final Euler error:   {error_euler[-1]:.10e}")
print(f"Final RK4 error:     {error_rk4[-1]:.10e}")


# Plot 1: Numerical and exact solutions
plt.figure(figsize=(9, 5))

plt.plot(x, y_exact, "k-", label="Exact solution", linewidth=2)
plt.plot(x, y_euler, "o--", label="Euler method")
plt.plot(x, y_rk4, "s:", label="RK4 method", linewidth=2)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Comparison of Euler and RK4 Methods")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("euler_vs_rk4_comparison.png", dpi=300)
plt.show()


# Plot 2: Absolute error comparison
plt.figure(figsize=(9, 5))

plt.semilogy(x, error_euler, "o-", label="Euler error")
plt.semilogy(x, error_rk4, "s-", label="RK4 error")

plt.xlabel("x")
plt.ylabel("Absolute error (log scale)")
plt.title("Absolute Error: Euler Method vs. RK4")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("euler_vs_rk4_comparison_error.png", dpi=300)
plt.show()

