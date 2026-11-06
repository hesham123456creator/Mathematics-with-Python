
import numpy as np
import matplotlib.pyplot as plt

# Differential equation
def f(x, y):
    return y - x**2 + 1

# Exact solution
def exact(x):
    return (x + 1)**2 - 0.5*np.exp(x)

# Initial conditions
x0 = 0
y0 = 0.5
h = 0.2
n = 10

x = np.linspace(x0, 2, n + 1)

# Euler's method
y_euler = np.zeros(n + 1)
y_euler[0] = y0

for i in range(n):
    y_euler[i+1] = (
        y_euler[i] + h*f(x[i], y_euler[i])
    )

# Fourth-order Runge-Kutta method
y_rk4 = np.zeros(n + 1)
y_rk4[0] = y0

for i in range(n):
    k1 = f(x[i], y_rk4[i])
    k2 = f(x[i] + h/2, y_rk4[i] + h*k1/2)
    k3 = f(x[i] + h/2, y_rk4[i] + h*k2/2)
    k4 = f(x[i] + h, y_rk4[i] + h*k3)

    y_rk4[i+1] = y_rk4[i] + (
        h/6 * (k1 + 2*k2 + 2*k3 + k4)
    )

# Exact values
y_exact = exact(x)

# Display results
print("x       Euler       RK4       Exact")

for i in range(n + 1):
    print(
        f"{x[i]:.1f}   "
        f"{y_euler[i]:.6f}   "
        f"{y_rk4[i]:.6f}   "
        f"{y_exact[i]:.6f}"
    )

# Calculate errors
error_euler = np.abs(y_euler - y_exact)
error_rk4 = np.abs(y_rk4 - y_exact)

print("\nMaximum Euler error:", np.max(error_euler))
print("Maximum RK4 error:", np.max(error_rk4))

# Plot solutions
plt.figure(figsize=(9, 5))

plt.plot(x, y_euler, "o-", label="Euler")
plt.plot(x, y_rk4, "s-", label="RK4")
plt.plot(x, y_exact, "--", label="Exact")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Euler vs RK4 Methods")

plt.legend()
plt.grid(True)
plt.show()
