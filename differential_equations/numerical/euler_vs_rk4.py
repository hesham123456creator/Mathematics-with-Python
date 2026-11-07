
import numpy as np
import matplotlib.pyplot as plt

# Differential equation: y' = -2y
def f(x, y):
    return -2 * y


# Exact solution
def exact_solution(x):
    return np.exp(-2 * x)


# Initial-value problem
x0 = 0.0
y0 = 1.0
x_end = 5.0
h = 0.1

# Discretize the interval
n_steps = int(round((x_end - x0) / h))
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


# Exact values
y_exact = exact_solution(x)


# Error analysis
error_euler = np.abs(y_euler - y_exact)
error_rk4 = np.abs(y_rk4 - y_exact)

print("Numerical solution of y' = -2y")
print("--------------------------------")
print(f"Step size: {h}")
print(f"Number of steps: {n_steps}")
print(f"Final x value: {x[-1]:.2f}")
print()
print("Maximum absolute errors:")
print(f"Euler method: {np.max(error_euler):.10e}")
print(f"RK4 method:   {np.max(error_rk4):.10e}")
print()
print("Final numerical values:")
print(f"Exact solution: {y_exact[-1]:.10e}")
print(f"Euler method:   {y_euler[-1]:.10e}")
print(f"RK4 method:     {y_rk4[-1]:.10e}")


# Plot 1: Comparison of numerical and exact solutions
plt.figure(figsize=(9, 5))

plt.plot(x, y_exact, label="Exact solution", linewidth=2)
plt.plot(x, y_euler, "--", label="Euler method")
plt.plot(x, y_rk4, ":", label="RK4 method", linewidth=2)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Euler Method vs. Fourth-Order Runge-Kutta")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("euler_vs_rk4.png", dpi=300)
plt.show()


# Plot 2: Absolute errors
plt.figure(figsize=(9, 5))

plt.semilogy(x, error_euler, label="Euler error")
plt.semilogy(x, error_rk4, label="RK4 error")

plt.xlabel("x")
plt.ylabel("Absolute error (log scale)")
plt.title("Error Comparison: Euler Method vs. RK4")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("euler_vs_rk4_error.png", dpi=300)
plt.show()

