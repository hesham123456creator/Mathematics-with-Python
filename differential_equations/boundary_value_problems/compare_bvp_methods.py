
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.linalg import solve_banded

# Boundary value problem:
# y'' + y = 0
# y(0) = 0, y(pi/2) = 1
# Exact solution: y(x) = sin(x)

# --------------------------------------------------
# 1. Finite Difference Method
# --------------------------------------------------

N = 40
a = 0.0
b = np.pi / 2

x_fd = np.linspace(a, b, N + 1)
h = (b - a) / N

# Construct the tridiagonal system
diagonal = (-2 + h**2) * np.ones(N - 1)
off_diagonal = np.ones(N - 2)

ab = np.zeros((3, N - 1))
ab[0, 1:] = off_diagonal
ab[1, :] = diagonal
ab[2, :-1] = off_diagonal

rhs = np.zeros(N - 1)

# Apply boundary conditions
rhs[-1] = -1.0

# Solve the linear system
y_interior = solve_banded((1, 1), ab, rhs)

y_fd = np.zeros(N + 1)
y_fd[1:N] = y_interior
y_fd[-1] = 1.0

# --------------------------------------------------
# 2. SciPy solve_bvp Method
# --------------------------------------------------

def ode(x, y):
    return np.vstack((y[1], -y[0]))


def bc(ya, yb):
    return np.array([ya[0], yb[0] - 1.0])


x_guess = np.linspace(a, b, 20)
y_guess = np.zeros((2, x_guess.size))

solution = solve_bvp(
    ode,
    bc,
    x_guess,
    y_guess,
    tol=1e-8
)

if not solution.success:
    raise RuntimeError(
        "solve_bvp failed: " + solution.message
    )

y_bvp = solution.sol(x_fd)[0]

# --------------------------------------------------
# 3. Exact Solution
# --------------------------------------------------

y_exact = np.sin(x_fd)

# --------------------------------------------------
# 4. Error Analysis
# --------------------------------------------------

error_fd = np.max(np.abs(y_fd - y_exact))
error_bvp = np.max(np.abs(y_bvp - y_exact))

print("Comparison of BVP Numerical Methods")
print("=" * 55)

print(f"Number of subintervals: {N}")
print(f"Step size: {h:.8f}")

print("\nFinite Difference Method")
print(f"Maximum absolute error: {error_fd:.10e}")

print("\nSciPy solve_bvp Method")
print(f"Maximum absolute error: {error_bvp:.10e}")
print(f"Number of mesh points: {solution.x.size}")
print(f"Solver successful: {solution.success}")

# --------------------------------------------------
# 5. Plot Comparison
# --------------------------------------------------

plt.figure(figsize=(9, 6))

plt.plot(
    x_fd,
    y_exact,
    label="Exact solution",
    linewidth=2.5
)

plt.plot(
    x_fd,
    y_fd,
    "o",
    label="Finite Difference Method",
    markersize=4
)

plt.plot(
    x_fd,
    y_bvp,
    "--",
    label="SciPy solve_bvp",
    linewidth=2
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Comparison of Numerical BVP Methods")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("compare_bvp_methods.png", dpi=300)

plt.show()
plt.close()

