
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Define the system of differential equations
def system(x, y):
    return np.vstack((
        y[1],
        -np.pi**2 * np.sin(np.pi * x)
    ))

# Define boundary conditions
def boundary_conditions(ya, yb):
    return np.array([
        ya[0],
        yb[0]
    ])

# Initial mesh
x = np.linspace(0, 1, 10)

# Initial guess for the solution
y_guess = np.zeros((2, x.size))

# Solve the boundary value problem
solution = solve_bvp(
    system,
    boundary_conditions,
    x,
    y_guess
)

# Check solver status
if not solution.success:
    raise RuntimeError(solution.message)

# Evaluate numerical solution
x_plot = np.linspace(0, 1, 200)
y_numerical = solution.sol(x_plot)[0]

# Exact solution
y_exact = np.sin(np.pi * x_plot)

# Calculate maximum error
error = np.max(np.abs(y_numerical - y_exact))

print("SciPy solve_bvp")
print("Solver successful:", solution.success)
print("Maximum absolute error:", error)

# Plot results
plt.figure(figsize=(8, 5))

plt.plot(
    x_plot,
    y_numerical,
    label="Numerical solution"
)

plt.plot(
    x_plot,
    y_exact,
    "--",
    label="Exact solution"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Boundary Value Problem Using solve_bvp")
plt.legend()
plt.grid(True)
plt.show()
