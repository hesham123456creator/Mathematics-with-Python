
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Nonlinear Bratu boundary-value problem:
# y'' + lambda * exp(y) = 0
# y(0) = 0, y(1) = 0

lam = 1.0

# Convert to a first-order system
def ode(x, y):
    return np.vstack((
        y[1],
        -lam * np.exp(y[0])
    ))

# Boundary conditions
def bc(ya, yb):
    return np.array([ya[0], yb[0]])

# Initial mesh
x = np.linspace(0, 1, 30)

# Nonzero initial guess
y_guess = np.zeros((2, x.size))
y_guess[0] = 0.1 * np.sin(np.pi * x)

# Solve the nonlinear BVP
solution = solve_bvp(ode, bc, x, y_guess)

if not solution.success:
    raise RuntimeError(solution.message)

# Evaluate the solution
x_plot = np.linspace(0, 1, 200)
y_plot = solution.sol(x_plot)[0]

print("Solver successful:", solution.success)
print("Solver message:", solution.message)
print("Maximum solution value:", np.max(y_plot))
print("Boundary value at x=0:", y_plot[0])
print("Boundary value at x=1:", y_plot[-1])

# Plot the nontrivial solution
plt.plot(x_plot, y_plot, label="Numerical solution")
plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonlinear Bratu Boundary Value Problem")
plt.legend()
plt.grid()
plt.savefig("nonlinear_bvp_nontrivial.png", dpi=300)
plt.show()
