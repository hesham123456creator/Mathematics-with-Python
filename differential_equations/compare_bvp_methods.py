import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# Problem: y'' + y = 0
# Boundary conditions: y(0) = 0, y(pi/2) = 1
# Exact solution: y(x) = sin(x)

# 1. Finite Difference Method
N = 40
x_fd = np.linspace(0, np.pi / 2, N + 1)
h = x_fd[1] - x_fd[0]

A = np.zeros((N - 1, N - 1))
rhs = np.zeros(N - 1)

for i in range(N - 1):
    A[i, i] = -2 + h**2

    if i > 0:
        A[i, i - 1] = 1

    if i < N - 2:
        A[i, i + 1] = 1

# Include the boundary condition y(pi/2) = 1
rhs[-1] = -1

y_fd = np.zeros(N + 1)
y_fd[1:N] = np.linalg.solve(A, rhs)
y_fd[-1] = 1

# 2. SciPy solve_bvp
def ode(x, y):
    return np.vstack((y[1], -y[0]))

def bc(ya, yb):
    return np.array([ya[0], yb[0] - 1])

x_guess = np.linspace(0, np.pi / 2, 20)
y_guess = np.zeros((2, x_guess.size))

solution = solve_bvp(ode, bc, x_guess, y_guess)

y_bvp = solution.sol(x_fd)[0]

# 3. Exact solution
y_exact = np.sin(x_fd)

# 4. Calculate errors
error_fd = np.max(np.abs(y_fd - y_exact))
error_bvp = np.max(np.abs(y_bvp - y_exact))

print("Finite Difference maximum error:", error_fd)
print("solve_bvp maximum error:", error_bvp)
print("solve_bvp successful:", solution.success)

# 5. Plot comparison
plt.plot(x_fd, y_exact, label="Exact solution")
plt.plot(x_fd, y_fd, "o", label="Finite Difference")
plt.plot(x_fd, y_bvp, "--", label="SciPy solve_bvp")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Comparison of BVP Numerical Methods")
plt.legend()
plt.grid()
plt.show()
