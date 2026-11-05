
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the differential equation
def f(x, y):
    return -2 * y

# Initial condition
y0 = [1]

# Interval
x_span = (0, 5)

# Points for numerical solution
x_eval = np.linspace(0, 5, 100)

# Solve the differential equation
solution = solve_ivp(
    f, x_span, y0, t_eval=x_eval
)

# Exact solution
y_exact = np.exp(-2 * x_eval)

# Display numerical results
print("   x       Numerical       Exact")
for i in range(0, 100, 10):
    print(
        f"{x_eval[i]:.2f}   "
        f"{solution.y[0][i]:.6f}   "
        f"{y_exact[i]:.6f}"
    )

# Plot both solutions
plt.plot(x_eval, solution.y[0],
         label="Numerical solution")

plt.plot(x_eval, y_exact, "--",
         label="Exact solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Numerical Solution of dy/dx = -2y")
plt.legend()
plt.grid(True)
plt.show()
