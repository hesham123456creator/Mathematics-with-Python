
import numpy as np
import matplotlib.pyplot as plt

# Exact solution
def exact_solution(x):
    return np.sin(np.pi * x)

# Solve using the finite difference method
def finite_difference(n):

    h = 1 / n
    x = np.linspace(0, 1, n + 1)

    # Number of interior points
    m = n - 1

    # Coefficient matrix
    A = np.zeros((m, m))

    for i in range(m):
        A[i, i] = -2

        if i > 0:
            A[i, i - 1] = 1

        if i < m - 1:
            A[i, i + 1] = 1

    A = A / h**2

    # Right-hand side
    b = -np.pi**2 * np.sin(np.pi * x[1:-1])

    # Solve the linear system
    y_interior = np.linalg.solve(A, b)

    # Include boundary conditions
    y = np.zeros(n + 1)
    y[1:-1] = y_interior

    # Calculate maximum error
    error = np.max(np.abs(y - exact_solution(x)))

    return x, y, error


# Different numbers of subintervals
n_values = [10, 20, 40, 80, 160]

errors = []

print("Accuracy comparison")
print("------------------------------")
print("n       Step size       Maximum error")

for n in n_values:

    x, y, error = finite_difference(n)

    h = 1 / n
    errors.append(error)

    print(f"{n:<7} {h:<15.6f} {error:.10e}")


# Plot error against step size
plt.figure(figsize=(8, 5))

h_values = [1 / n for n in n_values]

plt.loglog(h_values, errors, "o-", label="Numerical error")
plt.loglog(
    h_values,
    [errors[0] * (h / h_values[0])**2 for h in h_values],
    "--",
    label="Second-order reference"
)

plt.xlabel("Step size (h)")
plt.ylabel("Maximum absolute error")
plt.title("Finite Difference Method: Accuracy")
plt.legend()
plt.grid(True)
plt.show()
