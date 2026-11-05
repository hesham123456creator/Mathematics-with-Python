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
initial_conditions = [1, 0]

# Solve the system numerically
solution = solve_ivp(
    system,
    (0, 10),
    initial_conditions,
    method="RK45",
    dense_output=True
)

# Generate points for plotting
t = np.linspace(0, 10, 300)
x, y = solution.sol(t)

# Exact solutions
x_exact = np.cos(t)
y_exact = -np.sin(t)

# Calculate errors
error_x = np.max(np.abs(x - x_exact))
error_y = np.max(np.abs(y - y_exact))

print("Maximum error in x:", error_x)
print("Maximum error in y:", error_y)

# Plot both components
plt.plot(t, x, label="Numerical x(t)")
plt.plot(t, x_exact, "--", label="Exact x(t)")
plt.plot(t, y, label="Numerical y(t)")
plt.plot(t, y_exact, "--", label="Exact y(t)")

plt.xlabel("t")
plt.ylabel("Solutions")
plt.title("Numerical Solution of a System of ODEs")
plt.legend()
plt.grid()
plt.show()
