import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Nonlinear BVP:
# y'' + y^2 = 0
# y(0) = 0, y(1) = 0

# Convert to a first-order system
def ode(x, y):
    return np.vstack((y[1], -y[0]**2))

# Boundary conditions
def bc(ya, yb):
    return np.array([ya[0], yb[0]])

# Initial mesh
x = np.linspace(0, 1, 30)

# Initial guess
y_guess = np.zeros((2, x.size))

# Solve the nonlinear BVP
solution = solve_bvp(
    ode,
    bc,
    x,
    y_guess
)

# Evaluate the solution
x_plot = np.linspace(0, 1, 200)
y_plot = solution.sol(x_plot)[0]

print("Solver successful:", solution.success)
print("Solver message:", solution.message)
print("Maximum solution value:", np.max(np.abs(y_plot)))

# Plot solution
plt.plot(x_plot, y_plot, label="Numerical solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonlinear Boundary Value Problem")
plt.legend()
plt.grid()
plt.show()
