import numpy as np
import matplotlib.pyplot as plt

# Differential equation: y' = -2y
def f(x, y):
    return -2 * y

# Exact solution
def exact_solution(x):
    return np.exp(-2 * x)

# Initial conditions
x0 = 0
y0 = 1
h = 0.1
x_end = 5

# Discretize the interval
x = np.arange(x0, x_end + h, h)

# Euler's method
y_euler = np.zeros(len(x))
y_euler[0] = y0

for i in range(len(x) - 1):
    y_euler[i + 1] = y_euler[i] + h * f(x[i], y_euler[i])

# Fourth-order Runge-Kutta method
y_rk4 = np.zeros(len(x))
y_rk4[0] = y0

for i in range(len(x) - 1):
    k1 = f(x[i], y_rk4[i])
    k2 = f(x[i] + h/2, y_rk4[i] + h*k1/2)
    k3 = f(x[i] + h/2, y_rk4[i] + h*k2/2)
    k4 = f(x[i] + h, y_rk4[i] + h*k3)

    y_rk4[i + 1] = y_rk4[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

# Exact values
y_exact = exact_solution(x)

# Calculate maximum errors
error_euler = np.max(np.abs(y_euler - y_exact))
error_rk4 = np.max(np.abs(y_rk4 - y_exact))

print("Euler maximum error:", error_euler)
print("RK4 maximum error:", error_rk4)

# Plot comparison
plt.plot(x, y_exact, label="Exact solution")
plt.plot(x, y_euler, "--", label="Euler method")
plt.plot(x, y_rk4, ":", label="RK4 method")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Euler Method vs RK4")
plt.legend()
plt.grid()
plt.show()
