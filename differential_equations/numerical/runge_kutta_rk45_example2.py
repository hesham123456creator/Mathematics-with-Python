
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Differential equation: y' = y - x^2 + 1
def f(x, y):
    return y - x**2 + 1


# Exact analytical solution
def exact_solution(x):
    return (x + 1)**2 - 0.5 * np.exp(x)


# Initial-value problem
x_start = 0.0
x_end = 2.0
y0 = [0.5]

# Solver tolerances
relative_tolerance = 1e-8
absolute_tolerance = 1e-10

# Evaluation points
x_eval = np.linspace(x_start, x_end, 200)


# Numerical solution using adaptive RK45
solution = solve_ivp(
    f,
    (x_start, x_end),
    y0,
    method="RK45",
    t_eval=x_eval,
    rtol=relative_tolerance,
    atol=absolute_tolerance
)


# Validate solver output
if not solution.success:
    raise RuntimeError(
        f"RK45 solver failed: {solution.message}"
    )

if solution.y.shape[1] != len(x_eval):
    raise RuntimeError("The solver did not return all requested points.")


# Exact solution
y_exact = exact_solution(x_eval)
y_numerical = solution.y[0]

# Absolute errors
errors = np.abs(y_numerical - y_exact)


# Display selected numerical results
print("ADAPTIVE RK45 METHOD")
print("=" * 65)
print("Differential equation: y' = y - x^2 + 1")
print("Initial condition: y(0) = 0.5")
print(f"Interval: [{x_start}, {x_end}]")
print(f"Relative tolerance: {relative_tolerance:.1e}")
print(f"Absolute tolerance: {absolute_tolerance:.1e}")

print("\nNUMERICAL RESULTS")
print("-" * 65)
print(f"{'x':>8} {'RK45 solution':>18} "
      f"{'Exact solution':>18} {'Absolute error':>18}")

for i in range(0, len(x_eval), 20):
    print(
        f"{x_eval[i]:8.2f} "
        f"{y_numerical[i]:18.10f} "
        f"{y_exact[i]:18.10f} "
        f"{errors[i]:18.6e}"
    )

# Include the final point
i = len(x_eval) - 1
print(
    f"{x_eval[i]:8.2f} "
    f"{y_numerical[i]:18.10f} "
    f"{y_exact[i]:18.10f} "
    f"{errors[i]:18.6e}"
)


# Solver information
print("\nSOLVER PERFORMANCE")
print("-" * 40)
print("Solver success:", solution.success)
print("Accepted steps:", len(solution.t) - 1)
print("Function evaluations:", solution.nfev)
print("Number of Jacobian evaluations:", solution.njev)
print("Number of LU decompositions:", solution.nlu)

# Error summary
print("\nERROR ANALYSIS")
print("-" * 40)
print(f"Maximum absolute error: {np.max(errors):.12e}")
print(f"Final absolute error:   {errors[-1]:.12e}")


# Plot 1: Numerical and exact solutions
plt.figure(figsize=(9, 5))

plt.plot(
    x_eval, y_exact,
    "--", label="Exact solution", linewidth=2
)

plt.plot(
    x_eval, y_numerical,
    label="RK45 numerical solution"
)

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Adaptive RK45: Numerical and Exact Solutions")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("runge_kutta_rk45_example2.png", dpi=300)
plt.show()


# Plot 2: Absolute error
plt.figure(figsize=(9, 5))

plt.semilogy(
    x_eval, errors,
    label="RK45 absolute error"
)

plt.xlabel("x")
plt.ylabel("Absolute error (log scale)")
plt.title("RK45 Error Analysis")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("runge_kutta_rk45_example2_error.png", dpi=300)
plt.show()


# Plot 3: Adaptive step sizes
step_sizes = np.diff(solution.t)
step_locations = solution.t[:-1]

plt.figure(figsize=(9, 5))

plt.plot(step_locations, step_sizes, "o-", markersize=4)

plt.xlabel("x")
plt.ylabel("Accepted step size")
plt.title("Adaptive Step-Size Selection in RK45")
plt.grid(True)
plt.tight_layout()
plt.savefig("runge_kutta_rk45_example2_step_sizes.png", dpi=300)
plt.show()

