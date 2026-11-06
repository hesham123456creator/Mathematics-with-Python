
import numpy as np
import matplotlib.pyplot as plt

# Differential equation: y' = -2y
def ode(x, y):
    return -2 * y

# Initial condition
x0 = 0
y0 = 1

# Interval and step size
x_end = 5
h = 0.1
n = int((x_end - x0) / h)

x = np.linspace(x0, x_end, n + 1)
y = np.zeros(n + 1)
y[0] = y0

# Classical fourth-order Runge-Kutta method
for i in range(n):
    k1 = ode(x[i], y[i])
    k2 = ode(x[i] + h/2, y[i] + h*k1/2)
    k3 = ode(x[i] + h/2, y[i] + h*k2/2)
    k4 = ode(x[i] + h, y[i] + h*k3)

    y[i+1] = y[i] + (h/6) * (
        k1 + 2*k2 + 2*k3 + k4
    )

# Exact solution
y_exact = np.exp(-2*x)

# Calculate maximum absolute error
error = np.max(np.abs(y - y_exact))

print("Classical Runge-Kutta Fourth-Order Method")
print("Differential equation: y' = -2y")
print("Step size:", h)
print("Maximum absolute error:", error)

# Plot numerical and exact solutions
plt.figure(figsize=(9, 5))
plt.plot(x, y, "o-", label="RK4 numerical solution")
plt.plot(x, y_exact, "--", label="Exact solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Numerical Solution Using Classical RK4")
plt.legend()
plt.grid(True)
plt.show()
