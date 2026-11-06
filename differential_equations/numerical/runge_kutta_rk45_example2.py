
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the differential equation
def f(x, y):
    return y - x**2 + 1

# Initial condition
y0 = [0.5]

# Interval
x_span = (0, 2)

# Evaluation points
x_eval = np.linspace(0, 2, 100)

# Numerical solution using RK45
solution = solve_ivp(
    f,
    x_span,
    y0,
    method="RK45",
    t_eval=x_eval,
    rtol=1e-8,
    atol=1e-10
)

# Exact solution
y_exact = (
    (x_eval + 1)**2
    - 0.5*np.exp(x_eval)
)

# Display results
print("   x       Numerical       Exact")

for i in range(0, 100, 10):
    print(
        f"{x_eval[i]:.2f}   "
        f"{solution.y[0][i]:.8f}   "
        f"{y_exact[i]:.8f}"
    )

# Calculate maximum absolute error
error = np.max(
    np.abs(solution.y[0] - y_exact)
)

print("\nMaximum absolute error:", error)

# Plot both solutions
plt.figure(figsize=(9, 5))

plt.plot(
    x_eval, solution.y[0],
    label="RK45 Numerical Solution"
)

plt.plot(
    x_eval, y_exact, "--",
    label="Exact Solution"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Runge-Kutta RK45 Method")

plt.legend()
plt.grid(True)
plt.show()
