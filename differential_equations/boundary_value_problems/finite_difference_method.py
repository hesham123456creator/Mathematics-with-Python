import numpy as np
import matplotlib.pyplot as plt

# Boundary value problem:
# y''(x) = -pi^2 sin(pi*x)
# y(0) = 0, y(1) = 0

# Define the interval and number of subintervals
a = 0
b = 1
N = 20

h = (b - a) / N
x = np.linspace(a, b, N + 1)

# Construct the finite difference system
A = np.zeros((N - 1, N - 1))
rhs = np.zeros(N - 1)

for i in range(N - 1):
    A[i, i] = -2

    if i > 0:
        A[i, i - 1] = 1

    if i < N - 2:
        A[i, i + 1] = 1

    rhs[i] = -h**2 * np.pi**2 * np.sin(np.pi * x[i + 1])

# Solve the linear system
y_interior = np.linalg.solve(A, rhs)

# Include boundary values
y_numerical = np.zeros(N + 1)
y_numerical[1:N] = y_interior

# Exact solution
y_exact = np.sin(np.pi * x)

# Calculate maximum error
error = np.max(np.abs(y_numerical - y_exact))

print("Maximum absolute error:", error)

# Plot the results
plt.plot(x, y_numerical, "o-", label="Finite Difference")
plt.plot(x, y_exact, "--", label="Exact Solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Finite Difference Method for a BVP")
plt.legend()
plt.grid()
plt.show()
