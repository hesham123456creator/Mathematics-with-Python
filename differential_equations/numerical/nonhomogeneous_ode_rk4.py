
import numpy as np
import matplotlib.pyplot as plt

# Nonhomogeneous second-order differential equation:
# y'' + 3y' + 2y = exp(-x)
#
# Initial conditions:
# y(0) = 0, y'(0) = 0

# Convert the second-order equation into a system:
# y' = z
# z' = exp(-x) - 3z - 2y

def system(x, Y):
    y, z = Y

    dydx = z
    dzdx = np.exp(-x) - 3*z - 2*y

    return np.array([dydx, dzdx])


# Classical fourth-order Runge-Kutta method
def rk4(f, x0, Y0, x_end, h):
    n = int(round((x_end - x0) / h))

    x = np.linspace(x0, x_end, n + 1)
    Y = np.zeros((n + 1, len(Y0)))

    Y[0] = Y0

    for i in range(n):
        k1 = f(x[i], Y[i])
        k2 = f(x[i] + h/2, Y[i] + h*k1/2)
        k3 = f(x[i] + h/2, Y[i] + h*k2/2)
        k4 = f(x[i] + h, Y[i] + h*k3)

        Y[i+1] = Y[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

    return x, Y


# Initial conditions
x0 = 0
Y0 = [0, 0]

# Interval and step size
x_end = 5
h = 0.01

# Numerical solution
x, Y = rk4(system, x0, Y0, x_end, h)

y_numerical = Y[:, 0]

# Exact analytical solution
y_exact = (
    (x - 1)*np.exp(-x)
    + np.exp(-2*x)
)

# Calculate absolute errors
errors = np.abs(y_numerical - y_exact)
max_error = np.max(errors)

# Display results
print("Nonhomogeneous Second-Order ODE")
print("y'' + 3y' + 2y = exp(-x)")
print("Method: Classical RK4")
print("Step size:", h)

print("\n   x       Numerical       Exact       Absolute Error")

for i in range(0, len(x), 50):
    print(
        f"{x[i]:.2f}   "
        f"{y_numerical[i]:.8f}   "
        f"{y_exact[i]:.8f}   "
        f"{errors[i]:.2e}"
    )

print("\nMaximum absolute error:", max_error)

# Plot numerical and exact solutions
plt.figure(figsize=(9, 5))

plt.plot(
    x,
    y_numerical,
    label="Numerical solution (RK4)"
)

plt.plot(
    x,
    y_exact,
    "--",
    label="Exact solution"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Nonhomogeneous Second-Order ODE: RK4 Method")

plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
