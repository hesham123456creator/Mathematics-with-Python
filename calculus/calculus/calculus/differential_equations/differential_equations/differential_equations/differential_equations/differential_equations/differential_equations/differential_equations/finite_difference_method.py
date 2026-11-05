
import numpy as np
import matplotlib.pyplot as plt

# Number of subintervals
n = 20

# Step size
h = 1 / n

# Grid points
x = np.linspace(0, 1, n + 1)

# Interior grid points
x_interior = x[1:-1]

# Construct the coefficient matrix
A = np.zeros((n - 1, n - 1))

for i in range(n - 1):
    A[i, i] = -2

    if i > 0:
        A[i, i - 1] = 1

    if i < n - 2:
        A[i, i + 1] = 1

A = A / h**2

# Right-hand side
b = -np.pi**2 * np.sin(np.pi * x_interior)

# Solve the linear system
y_interior = np.linalg.solve(A, b)

# Include boundary conditions
y_numerical = np.zeros(n + 1)
y_numerical[1:-1] = y_interior

# Exact solution
y_exact = np.sin(np.pi * x)

# Calculate maximum absolute error
error = np.max(np.abs(y_numerical - y_exact))

print("Finite Difference Method")
print("Maximum absolute error:", error)

# Plot numerical and exact solutions
plt.figure(figsize=(8, 5))

plt.plot(x, y_numerical, "o-", label="Numerical solution")
plt.plot(x, y_exact, "--", label="Exact solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Finite Difference Method")
plt.legend()
plt.grid(True)
plt.show()
