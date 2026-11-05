
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Nonhomogeneous second-order ODE:
# y'' + 3y' + 2y = exp(-x)

def system(x, Y):
    y, z = Y
    return [z, np.exp(-x) - 3*z - 2*y]

# Initial conditions
y0 = [0, 0]

# Interval
x_span = (0, 5)

# Evaluation points
x_eval = np.linspace(0, 5, 100)

# Numerical solution
solution = solve_ivp(
    system, x_span, y0, t_eval=x_eval
)

# Exact solution
y_exact = (
    np.exp(-x_eval)
    - np.exp(-2*x_eval)
    + x_eval*np.exp(-x_eval)
)

# Display results
print("   x       Numerical       Exact")

for i in range(0, 100, 10):
    print(
        f"{x_eval[i]:.2f}   "
        f"{solution.y[0][i]:.6f}   "
        f"{y_exact[i]:.6f}"
    )

# Plot solutions
plt.figure(figsize=(9, 5))

plt.plot(x_eval, solution.y[0],
         label="Numerical solution")

plt.plot(x_eval, y_exact, "--",
         label="Exact solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonhomogeneous Second-Order ODE")
plt.legend()
plt.grid(True)
plt.show()
