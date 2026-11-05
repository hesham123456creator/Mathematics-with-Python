
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Define the nonlinear system
def system(x, y, p):
    lam = p[0]

    return np.vstack((
        y[1],
        -lam * y[0]**2
    ))

# Define boundary conditions
def boundary_conditions(ya, yb, p):
    return np.array([
        ya[0],       # y(0) = 0
        ya[1] - 1,   # y'(0) = 1
        yb[0] - 1    # y(1) = 1
    ])

# Initial mesh
x = np.linspace(0, 1, 20)

# Initial guess for the solution
y_guess = np.zeros((2, x.size))
y_guess[0] = x
y_guess[1] = 1

# Initial guess for the unknown parameter
p_guess = np.array([1.0])

# Solve the boundary value problem
solution = solve_bvp(
    system,
    boundary_conditions,
    x,
    y_guess,
    p=p_guess,
    tol=1e-6,
    max_nodes=1000
)

# Check convergence
if not solution.success:
    raise RuntimeError(solution.message)

# Extract the parameter
lam = solution.p[0]

# Evaluate the solution
x_plot = np.linspace(0, 1, 200)
y_numerical = solution.sol(x_plot)[0]

# Display results
print("Nonlinear BVP with an unknown parameter")
print("---------------------------------------")
print("Solver successful:", solution.success)
print("Estimated parameter lambda:", lam)
print("y(0):", solution.sol(0)[0])
print("y'(0):", solution.sol(0)[1])
print("y(1):", solution.sol(1)[0])

# Plot numerical solution
plt.figure(figsize=(8, 5))

plt.plot(
    x_plot,
    y_numerical,
    label="Numerical solution"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonlinear BVP with Unknown Parameter")
plt.legend()
plt.grid(True)
plt.show()
