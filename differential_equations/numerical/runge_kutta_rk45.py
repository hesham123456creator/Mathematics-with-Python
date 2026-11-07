
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Differential equation: y' = -2y
def ode(x, y):
    return -2 * y


# Exact analytical solution
def exact_solution(x):
    return np.exp(-2 * x)


# Initial-value problem
x_start = 0.0
x_end = 5.0
y0 = [1.0]

# RK45 accuracy settings
relative_tolerance = 1e-8
absolute_tolerance = 1e-10


# Numerical solution using adaptive RK45
solution = solve_ivp(
    ode,
    (x_start, x_end),
    y0,
    method="RK45",
    rtol=relative_tolerance,
    atol=absolute_tolerance,
    dense_output=True
)


# Check solver status
if not solution.success:
    raise RuntimeError(
        f"RK45 solver failed: {solution.message}"
    )


# Generate points for evaluating the solution
x = np.linspace(x_start, x_end, 500)
y_numerical = solution.sol(x)[0]

# Exact solution
y_exact = exact_solution(x)

# Absolute errors
errors = np.abs(y_numerical - y_exact)


# Display numerical results
print("ADAPTIVE RUNGE-KUTTA-FEHLBERG METHOD (RK45)")
print("=" * 55)
print("Differential equation: y' = -2y")
print("Initial condition: y(0) = 1")
print(f"Interval: [{x_start}, {x_end}]")
print(f"Relative tolerance: {relative_tolerance:.1e}")
print(f"Absolute tolerance: {absolute_tolerance:.1e}")

print("\nSOLVER INFORMATION")
print("-" * 40)
print("Solver success:", solution.success)
print("Solver message:", solution.message)
print("Accepted steps:", len(solution.t) - 1)
print("Function evaluations:", solution.nfev)

print("\nERROR ANALYSIS")
print("-" * 40)
print(f"Maximum absolute error: {np.max(errors):.12e}")
print(f"Final absolute error:   {errors[-1]:.12e}")
print(f"Final numerical value:  {y_numerical[-1]:.12e}")
print(f"Final exact value:      {y_exact[-1]:.12e}")


# Plot 1: Numerical and exact solutions
plt.figure(figsize=(9, 5))

plt.plot(x, y_exact, "--", label="Exact solution", linewidth=2)
plt.plot(x, y_numerical, label="RK45 numerical solution")

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Numerical Solution Using Adaptive RK45")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("runge_kutta_rk45.png", dpi=300)
plt.show()


# Plot 2: Absolute error
plt.figure(figsize=(9, 5))

plt.semilogy(x, errors, label="RK45 absolute error")

plt.xlabel("x")
plt.ylabel("Absolute error (log scale)")
plt.title("RK45 Numerical Error")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("runge_kutta_rk45_error.png", dpi=300)
plt.show()


# Plot 3: Adaptive step sizes
step_sizes = np.diff(solution.t)
step_locations = solution.t[:-1]

plt.figure(figsize=(9, 5))

plt.plot(step_locations, step_sizes, "o-", markersize=3)

plt.xlabel("x")
plt.ylabel("Accepted step size")
plt.title("Adaptive Step-Size Selection in RK45")
plt.grid(True)
plt.tight_layout()
plt.savefig("runge_kutta_rk45_step_sizes.png", dpi=300)
plt.show()
