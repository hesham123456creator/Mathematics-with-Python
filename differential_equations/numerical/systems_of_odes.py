

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# System of differential equations:
# x' = y
# y' = -x

def system(t, z):
    x, y = z
    return [y, -x]


# Exact solutions
def exact_x(t):
    return np.cos(t)


def exact_y(t):
    return -np.sin(t)


# Initial conditions and interval
t_start = 0.0
t_end = 10.0
initial_conditions = [1.0, 0.0]

# Evaluation points
t_eval = np.linspace(t_start, t_end, 300)


# Numerical solution using adaptive RK45
solution = solve_ivp(
    system,
    (t_start, t_end),
    initial_conditions,
    method="RK45",
    t_eval=t_eval,
    dense_output=True,
    rtol=1e-9,
    atol=1e-11
)


# Validate numerical solution
if not solution.success:
    raise RuntimeError(
        f"Solver failed: {solution.message}"
    )

# Numerical components
x_numerical = solution.y[0]
y_numerical = solution.y[1]

# Exact components
x_true = exact_x(t_eval)
y_true = exact_y(t_eval)


# Absolute errors
error_x = np.abs(x_numerical - x_true)
error_y = np.abs(y_numerical - y_true)


# Energy conservation:
# E = (x^2 + y^2)/2
energy = 0.5 * (x_numerical**2 + y_numerical**2)
exact_energy = 0.5 * (1.0**2 + 0.0**2)
energy_error = np.abs(energy - exact_energy)


# Numerical results
print("SYSTEM OF DIFFERENTIAL EQUATIONS")
print("=" * 65)
print("x' = y")
print("y' = -x")
print("Initial conditions: x(0) = 1, y(0) = 0")
print("Numerical method: Adaptive RK45")

print("\nSOLVER PERFORMANCE")
print("-" * 40)
print("Solver success:", solution.success)
print("Solver message:", solution.message)
print("Accepted steps:", len(solution.t) - 1)
print("Function evaluations:", solution.nfev)

print("\nERROR ANALYSIS")
print("-" * 40)
print(f"Maximum error in x: {np.max(error_x):.12e}")
print(f"Maximum error in y: {np.max(error_y):.12e}")
print(f"Maximum energy error: {np.max(energy_error):.12e}")

print("\nSELECTED NUMERICAL VALUES")
print("-" * 65)
print(f"{'t':>8} {'Numerical x':>16} {'Exact x':>16} "
      f"{'Numerical y':>16}")

for i in range(0, len(t_eval), 30):
    print(
        f"{t_eval[i]:8.2f} "
        f"{x_numerical[i]:16.10f} "
        f"{x_true[i]:16.10f} "
        f"{y_numerical[i]:16.10f}"
    )


# Plot 1: Numerical and exact solutions
plt.figure(figsize=(10, 5))

plt.plot(t_eval, x_numerical, label="Numerical x(t)")
plt.plot(t_eval, x_true, "--", label="Exact x(t)")

plt.plot(t_eval, y_numerical, label="Numerical y(t)")
plt.plot(t_eval, y_true, "--", label="Exact y(t)")

plt.xlabel("t")
plt.ylabel("Solutions")
plt.title("Numerical Solution of a System of ODEs")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("systems_of_odes.png", dpi=300)
plt.show()


# Plot 2: Absolute errors
plt.figure(figsize=(9, 5))

plt.semilogy(t_eval, error_x, label="Error in x(t)")
plt.semilogy(t_eval, error_y, label="Error in y(t)")

plt.xlabel("t")
plt.ylabel("Absolute error (log scale)")
plt.title("Numerical Error Analysis")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("systems_of_odes_error.png", dpi=300)
plt.show()


# Plot 3: Phase-plane trajectory
plt.figure(figsize=(7, 7))

plt.plot(x_numerical, y_numerical, label="Numerical trajectory")
plt.plot(x_true, y_true, "--", label="Exact trajectory")

plt.scatter(
    [initial_conditions[0]],
    [initial_conditions[1]],
    label="Initial point"
)

plt.xlabel("x(t)")
plt.ylabel("y(t)")
plt.title("Phase-Plane Trajectory")
plt.axis("equal")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("systems_of_odes_phase_plane.png", dpi=300)
plt.show()


# Plot 4: Energy conservation
plt.figure(figsize=(9, 5))

plt.plot(t_eval, energy, label="Numerical energy")
plt.axhline(
    exact_energy,
    linestyle="--",
    label="Exact constant energy"
)

plt.xlabel("t")
plt.ylabel("Energy")
plt.title("Energy Conservation in the ODE System")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("systems_of_odes_energy.png", dpi=300)
plt.show()
