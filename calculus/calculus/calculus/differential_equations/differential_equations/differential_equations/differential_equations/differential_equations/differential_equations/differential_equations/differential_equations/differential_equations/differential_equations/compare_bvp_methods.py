
import numpy as np
import matplotlib.pyplot as plt
import time

from scipy.integrate import solve_bvp

# Exact solution
def exact_solution(x):
    return np.sin(np.pi * x)


# -----------------------------------
# Method 1: Finite Difference Method
# -----------------------------------

def finite_difference(n):

    h = 1 / n
    x = np.linspace(0, 1, n + 1)

    m = n - 1
    A = np.zeros((m, m))

    for i in range(m):
        A[i, i] = -2

        if i > 0:
            A[i, i - 1] = 1

        if i < m - 1:
            A[i, i + 1] = 1

    A = A / h**2

    b = -np.pi**2 * np.sin(np.pi * x[1:-1])

    y_interior = np.linalg.solve(A, b)

    y = np.zeros(n + 1)
    y[1:-1] = y_interior

    return x, y


# -----------------------------------
# Method 2: SciPy solve_bvp
# -----------------------------------

def scipy_bvp(n):

    def system(x, y):
        return np.vstack((
            y[1],
            -np.pi**2 * np.sin(np.pi * x)
        ))

    def boundary_conditions(ya, yb):
        return np.array([ya[0], yb[0]])

    x = np.linspace(0, 1, n + 1)
    y_guess = np.zeros((2, x.size))

    solution = solve_bvp(
        system,
        boundary_conditions,
        x,
        y_guess,
        tol=1e-6
    )

    if not solution.success:
        raise RuntimeError(solution.message)

    return solution


# -----------------------------------
# Compare both methods
# -----------------------------------

n = 100

# Finite Difference Method
start = time.perf_counter()
x_fd, y_fd = finite_difference(n)
fd_time = time.perf_counter() - start

fd_error = np.max(
    np.abs(y_fd - exact_solution(x_fd))
)

# SciPy solve_bvp
start = time.perf_counter()
bvp_solution = scipy_bvp(n)
bvp_time = time.perf_counter() - start

x_plot = np.linspace(0, 1, 1000)
y_bvp = bvp_solution.sol(x_plot)[0]

bvp_error = np.max(
    np.abs(y_bvp - exact_solution(x_plot))
)


# -----------------------------------
# Display results
# -----------------------------------

print("Comparison of Numerical Methods")
print("----------------------------------------")

print("Finite Difference Method")
print("Maximum error:", fd_error)
print("Computational time:", fd_time, "seconds")

print("\nSciPy solve_bvp")
print("Maximum error:", bvp_error)
print("Computational time:", bvp_time, "seconds")


# -----------------------------------
# Plot numerical solutions
# -----------------------------------

plt.figure(figsize=(9, 5))

plt.plot(
    x_fd,
    y_fd,
    "o",
    markersize=3,
    label="Finite Difference"
)

plt.plot(
    x_plot,
    y_bvp,
    label="SciPy solve_bvp"
)

plt.plot(
    x_plot,
    exact_solution(x_plot),
    "--",
    label="Exact solution"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Comparison of BVP Numerical Methods")
plt.legend()
plt.grid(True)
plt.show()


# -----------------------------------
# Plot absolute errors
# -----------------------------------

plt.figure(figsize=(9, 5))

plt.semilogy(
    x_fd,
    np.abs(y_fd - exact_solution(x_fd)),
    "o-",
    markersize=3,
    label="Finite Difference Error"
)

plt.semilogy(
    x_plot,
    np.abs(y_bvp - exact_solution(x_plot)),
    label="SciPy solve_bvp Error"
)

plt.xlabel("x")
plt.ylabel("Absolute Error")
plt.title("Error Comparison")
plt.legend()
plt.grid(True)
plt.show()
