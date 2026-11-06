import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Convert the second-order equation into a system:
# y1' = y2
# y2' = -y1

def ode(x, y):
    return np.vstack((y[1], -y[0]))

# Boundary conditions
def boundary_conditions(ya, yb):
    return np.array([ya[0], yb[0] - 1])

# Initial mesh and guess
x = np.linspace(0, np.pi / 2, 20)
y_guess = np.zeros((2, x.size))

# Solve the boundary value problem
solution = solve_bvp(
    ode,
    boundary_conditions,
    x,
    y_guess
)

# Evaluate the numerical solution
x_plot = np.linspace(0, np.pi / 2, 200)
y_numerical = solution.sol(x_plot)[0]

# Exact solution
y_exact = np.sin(x_plot)

# Calculate maximum error
error = np.max(np.abs(y_numerical - y_exact))

print("Solver successful:", solution.success)
print("Maximum absolute error:", error)

# Plot results
plt.plot(x_plot, y_numerical, label="Numerical solution")
plt.plot(x_plot, y_exact, "--", label="Exact solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Boundary Value Problem Using SciPy")
plt.legend()
plt.grid()
plt.show()
