
import numpy as np
import matplotlib.pyplot as plt

# First-order differential equation:
# dy/dx = -2y
# Initial condition: y(0) = 1

def f(x, y):
    return -2 * y


def rk4(f, x0, y0, x_end, h):
    n = int(round((x_end - x0) / h))

    x = np.linspace(x0, x_end, n + 1)
    y = np.zeros(n + 1)
    y[0] = y0

    for i in range(n):
        k1 = f(x[i], y[i])
        k2 = f(x[i] + h/2, y[i] + h*k1/2)
        k3 = f(x[i] + h/2, y[i] + h*k2/2)
        k4 = f(x[i] + h, y[i] + h*k3)

        y[i+1] = y[i] + (h/6) * (
            k1 + 2*k2 + 2*k3 + k4
        )

    return x, y


# Initial condition and numerical parameters
x0 = 0
y0 = 1
x_end = 5
h = 0.01

# Compute numerical solution using RK4
x, y_numerical = rk4(f, x0, y0, x_end, h)

# Exact solution
y_exact = np.exp(-2*x)

# Calculate absolute errors
errors = np.abs(y_numerical - y_exact)
max_error = np.max(errors)

# Display numerical results
print("First-Order Differential Equation")
print("dy/dx = -2y, y(0) = 1")
print("Method: Classical RK4")
print(f"Step size: {h}")
print("\n   x       Numerical       Exact       Absolute Error")

for i in range(0, len(x), 50):
    print(
        f"{x[i]:.2f}   "
        f"{y_numerical[i]:.8f}   "
        f"{y_exact[i]:.8f}   "
        f"{errors[i]:.2e}"
    )

print(f"\nMaximum absolute error: {max_error:.12e}")

# Plot numerical and exact solutions
plt.plot(x, y_numerical, label="Numerical RK4 solution")
plt.plot(x, y_exact, "--", label="Exact solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("First-Order ODE: Classical RK4 Method")
plt.legend()
plt.grid(True)
plt.show()
