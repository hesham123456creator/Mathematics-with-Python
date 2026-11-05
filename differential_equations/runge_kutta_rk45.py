import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Differential equation: y' = -2y
def ode(x, y):
    return -2 * y

# Initial condition
x0 = 0
y0 = [1]

# Interval
x_start = 0
x_end = 5

# Numerical solution using RK45
solution = solve_ivp(
    ode,
    (x_start, x_end),
    y0,
    method="RK45",
    dense_output=True
)

# Generate points for plotting
x = np.linspace(x_start, x_end, 200)
y_numerical = solution.sol(x)[0]

# Exact solution
y_exact = np.exp(-2 * x)

# Calculate maximum absolute error
error = np.max(np.abs(y_numerical - y_exact))

print("Number of function evaluations:", solution.nfev)
print("Maximum absolute error:", error)

# Plot numerical and exact solutions
plt.plot(x, y_numerical, label="RK45 numerical solution")
plt.plot(x, y_exact, "--", label="Exact solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Numerical Solution Using RK45")
plt.legend()
plt.grid()
plt.show()
