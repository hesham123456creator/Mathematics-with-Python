
import numpy as np
import matplotlib.pyplot as plt

# System of differential equations:
# x' = y
# y' = -x

def system(t, z):
    x, y = z
    return np.array([y, -x])

# Initial conditions
t0 = 0
z0 = np.array([1.0, 0.0])

# Time interval and step size
t_end = 10
h = 0.01
n = int((t_end - t0) / h)

t = np.linspace(t0, t_end, n + 1)
z = np.zeros((n + 1, 2))
z[0] = z0

# Classical fourth-order Runge-Kutta method
for i in range(n):
    k1 = system(t[i], z[i])
    k2 = system(t[i] + h/2, z[i] + h*k1/2)
    k3 = system(t[i] + h/2, z[i] + h*k2/2)
    k4 = system(t[i] + h, z[i] + h*k3)

    z[i+1] = z[i] + (h/6) * (
        k1 + 2*k2 + 2*k3 + k4
    )

# Numerical solutions
x_numerical = z[:, 0]
y_numerical = z[:, 1]

# Exact solutions
x_exact = np.cos(t)
y_exact = -np.sin(t)

# Calculate maximum errors
error_x = np.max(np.abs(x_numerical - x_exact))
error_y = np.max(np.abs(y_numerical - y_exact))

print("System of Differential Equations")
print("x' = y, y' = -x")
print("Method: Classical RK4")
print("Step size:", h)
print("Maximum error in x:", error_x)
print("Maximum error in y:", error_y)

# Plot the solutions
plt.figure(figsize=(9, 5))

plt.plot(t, x_numerical, label="Numerical x(t)")
plt.plot(t, y_numerical, label="Numerical y(t)")
plt.plot(t, x_exact, "--", label="Exact x(t)")
plt.plot(t, y_exact, "--", label="Exact y(t)")

plt.xlabel("t")
plt.ylabel("Solution")
plt.title("System of ODEs Using Classical RK4")
plt.legend()
plt.grid(True)
plt.show()
