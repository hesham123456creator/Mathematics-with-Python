
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Define the nonlinear differential equation
def system(x, y):
    return np.vstack((
        y[1],
        -y[0]**2
    ))

# Define boundary conditions
def boundary_conditions(ya, yb):
    return np.array([
        ya[0],
        yb[0] - 1
    ])

# Initial mesh
x = np.linspace(0, 1, 20)

# Initial guess
y_guess = np.zeros((2, x.size))

# Improve the initial guess for y
y_guess[0] = x

# Solve the nonlinear boundary value problem
solution = solve_bvp(
    system,
    boundary_conditions,
    x,
    y_guess,
    tol=1e-6,
    max_nodes=1000
)

# Check convergence
if not solution.success:
    raise RuntimeError(solution.message)

# Evaluate numerical solution
x_plot = np.linspace(0, 1, 200)
y_numerical = solution.sol(x_plot)[0]

# Check boundary values
print("Nonlinear BVP")
print("Solver successful:", solution.success)
print("y(0) =", solution.sol(0)[0])
print("y(1) =", solution.sol(1)[0])
print("Number of mesh points:", solution.x.size)

# Plot solution
plt.figure(figsize=(8, 5))

plt.plot(
    x_plot,
    y_numerical,
    label="Numerical solution"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonlinear Boundary Value Problem")
plt.legend()
plt.grid(True)
plt.show()
