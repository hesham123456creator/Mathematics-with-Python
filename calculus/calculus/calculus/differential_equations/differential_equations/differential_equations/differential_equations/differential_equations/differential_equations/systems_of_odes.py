
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the system of differential equations
def system(t, z):
    x, y = z
    dxdt = y
    dydt = -x
    return [dxdt, dydt]

# Initial conditions
initial_conditions = [1, 0]

# Time interval
t_span = (0, 10)
t_values = np.linspace(0, 10, 200)

# Solve the system
solution = solve_ivp(
    system,
    t_span,
    initial_conditions,
    t_eval=t_values
)

# Numerical solutions
x_numerical = solution.y[0]
y_numerical = solution.y[1]

# Exact solutions
x_exact = np.cos(t_values)
y_exact = -np.sin(t_values)

# Plot the solutions
plt.figure(figsize=(9, 5))

plt.plot(t_values, x_numerical, label="Numerical x(t)")
plt.plot(t_values, y_numerical, label="Numerical y(t)")
plt.plot(t_values, x_exact, "--", label="Exact x(t)")
plt.plot(t_values, y_exact, "--", label="Exact y(t)")

plt.xlabel("t")
plt.ylabel("Solution")
plt.title("System of Differential Equations")
plt.legend()
plt.grid(True)
plt.show()
