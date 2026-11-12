
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
        ya[0],          # y(0) = 0
        ya[1] - 1,      # y'(0) = 1
        yb[0] - 0.5     # y(1) = 0.5
    ])


# Initial mesh
x = np.linspace(0, 1, 20)

# Initial guess for the solution
y_guess = np.zeros((2, x.size))
y_guess[0] = x - 0.5 * x**2
y_guess[1] = 1 - x

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

# Solver's collocation residual estimates
max_collocation_residual = np.max(solution.rms_residuals)

# Extract the estimated parameter
lam = solution.p[0]

# Evaluate the numerical solution
x_plot = np.linspace(0, 1, 200)
y_values = solution.sol(x_plot)
y_numerical = y_values[0]

# Evaluate the second derivative
y_second = solution.sol(x_plot, 2)[0]

# Calculate the ODE residual
ode_residual = y_second + lam * y_numerical**2
max_ode_residual = np.max(np.abs(ode_residual))

# Calculate boundary values
y0 = solution.sol(0)[0]
yp0 = solution.sol(0)[1]
y1 = solution.sol(1)[0]

# Calculate boundary residuals
left_residual = abs(y0)
initial_slope_residual = abs(yp0 - 1)
right_residual = abs(y1 - 0.5)

# Calculate maximum solution value
max_solution = np.max(y_numerical)


# Display results
print("Nonlinear BVP with an unknown parameter")
print("=" * 45)

print("\nSolver Information")
print("Solver successful:", solution.success)
print("Solver message:", solution.message)
print("Number of mesh points:", solution.x.size)

print("\nEstimated Parameter")
print(f"Lambda: {lam:.12f}")

print("\nBoundary Conditions")
print(f"y(0): {y0:.12f}")
print(f"y'(0): {yp0:.12f}")
print(f"y(1): {y1:.12f}")

print("\nAccuracy Analysis")
print(f"Maximum solution value: {max_solution:.12f}")
print(f"ODE residual: {max_ode_residual:.6e}")
print(
    f"Maximum collocation residual: "
    f"{max_collocation_residual:.6e}"
)

print("\nBoundary Residuals")
print(f"y(0) residual: {left_residual:.6e}")
print(f"y'(0) residual: {initial_slope_residual:.6e}")
print(f"y(1) residual: {right_residual:.6e}")


# Plot numerical solution
plt.figure(figsize=(8, 5))

plt.plot(
    x_plot,
    y_numerical,
    label="Numerical solution",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonlinear BVP with Unknown Parameter")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
