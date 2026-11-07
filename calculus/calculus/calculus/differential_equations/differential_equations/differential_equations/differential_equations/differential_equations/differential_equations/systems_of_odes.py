
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the system:
# x' = y
# y' = -x
def system(t, z):
    x, y = z
    return [y, -x]


# Initial conditions
initial_conditions = [1.0, 0.0]

# Time interval
t_span = (0.0, 10.0)
t_values = np.linspace(*t_span, 300)

# Solve using adaptive RK45
solution = solve_ivp(
    system,
    t_span,
    initial_conditions,
    t_eval=t_values,
    method="RK45",
    rtol=1e-9,
    atol=1e-11
)

# Check solver success
if not solution.success:
    raise RuntimeError(solution.message)

# Numerical solutions
x_numerical = solution.y[0]
y_numerical = solution.y[1]

# Exact solutions
x_exact = np.cos(t_values)
y_exact = -np.sin(t_values)

# Absolute errors
x_error = np.abs(x_numerical - x_exact)
y_error = np.abs(y_numerical - y_exact)

# Total energy: E = (x^2 + y^2) / 2
energy = 0.5 * (x_numerical**2 + y_numerical**2)

# Display numerical results
print("System of Differential Equations")
print("-" * 45)
print("Method: Adaptive Runge-Kutta (RK45)")
print("Solver message:", solution.message)
print("Accepted time steps:", len(solution.t))
print("Function evaluations:", solution.nfev)
print("Maximum error in x(t):", f"{np.max(x_error):.3e}")
print("Maximum error in y(t):", f"{np.max(y_error):.3e}")
print("Maximum energy deviation:", f"{np.max(np.abs(energy - 0.5)):.3e}")

# Plot numerical and exact solutions
plt.figure(figsize=(10, 6))
plt.plot(t_values, x_numerical, label="Numerical x(t)")
plt.plot(t_values, y_numerical, label="Numerical y(t)")
plt.plot(t_values, x_exact, "--", label="Exact x(t)")
plt.plot(t_values, y_exact, "--", label="Exact y(t)")
plt.xlabel("Time")
plt.ylabel("Solution")
plt.title("Numerical and Exact Solutions")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("systems_of_odes_comparison.png", dpi=300)
plt.show()

# Plot absolute errors
plt.figure(figsize=(10, 5))
plt.semilogy(t_values, x_error, label="Error in x(t)")
plt.semilogy(t_values, y_error, label="Error in y(t)")
plt.xlabel("Time")
plt.ylabel("Absolute error")
plt.title("Numerical Errors")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("systems_of_odes_error.png", dpi=300)
plt.show()

# Plot energy conservation
plt.figure(figsize=(10, 5))
plt.plot(t_values, energy, label="Numerical energy")
plt.axhline(0.5, linestyle="--", label="Exact energy")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("Energy Conservation")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("systems_of_odes_energy.png", dpi=300)
plt.show()
