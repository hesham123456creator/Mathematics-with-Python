
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the system of first-order equations
def system(x, Y):
    y1, y2 = Y
    return [y2, -3*y2 - 2*y1]

# Initial conditions
y0 = [1, 0]

# Interval
x_span = (0, 5)

# Evaluation points
x_eval = np.linspace(0, 5, 100)

# Numerical solution
solution = solve_ivp(
    system,
    x_span,
    y0,
    t_eval=x_eval
)

# Exact solution
y_exact = 2*np.exp(-x_eval) - np.exp(-2*x_eval)

# Display numerical results
print("   x       Numerical       Exact")

for i in range(0, 100, 10):
    print(
        f"{x_eval[i]:.2f}   "
        f"{solution.y[0][i]:.6f}   "
        f"{y_exact[i]:.6f}"
    )

# Plot numerical and exact solutions
plt.figure(figsize=(9, 5))

plt.plot(
    x_eval, solution.y[0],
    label="Numerical solution"
)

plt.plot(
    x_eval, y_exact, "--",
    label="Exact solution"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Second-Order Differential Equation")

plt.legend()
plt.grid(True)
plt.show()
