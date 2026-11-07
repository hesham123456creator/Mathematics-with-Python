
import numpy as np
import matplotlib.pyplot as plt


# System of differential equations:
# x' = y
# y' = -x

def system(t, z):
    x, y = z
    return np.array([y, -x])


# Exact solutions
def exact_x(t):
    return np.cos(t)


def exact_y(t):
    return -np.sin(t)


# Classical fourth-order Runge-Kutta method
def rk4(f, t0, z0, t_end, h):

    if h <= 0:
        raise ValueError("Step size must be positive.")

    if t_end <= t0:
        raise ValueError("Final time must exceed initial time.")

    n = int(np.ceil((t_end - t0) / h))
    t = np.linspace(t0, t_end, n + 1)
    h_actual = (t_end - t0) / n

    z = np.zeros((n + 1, len(z0)))
    z[0] = z0

    for i in range(n):

        hi = t[i + 1] - t[i]

        k1 = f(t[i], z[i])
        k2 = f(t[i] + hi/2, z[i] + hi*k1/2)
        k3 = f(t[i] + hi/2, z[i] + hi*k2/2)
        k4 = f(t[i] + hi, z[i] + hi*k3)

        z[i + 1] = z[i] + (hi/6) * (
            k1 + 2*k2 + 2*k3 + k4
        )

    return t, z, h_actual


# Initial conditions and interval
t0 = 0.0
t_end = 10.0
z0 = np.array([1.0, 0.0])
h = 0.01

# Numerical solution
t, z, h_actual = rk4(system, t0, z0, t_end, h)

x_numerical = z[:, 0]
y_numerical = z[:, 1]

# Exact solutions
x_true = exact_x(t)
y_true = exact_y(t)

# Absolute errors
errors_x = np.abs(x_numerical - x_true)
errors_y = np.abs(y_numerical - y_true)


# Energy conservation:
# E = (x^2 + y^2)/2
energy = 0.5 * (x_numerical**2 + y_numerical**2)
exact_energy = 0.5
energy_error = np.abs(energy - exact_energy)


# Display numerical results
print("SYSTEM OF DIFFERENTIAL EQUATIONS")
print("=" * 65)
print("x' = y")
print("y' = -x")
print("Initial conditions: x(0) = 1, y(0) = 0")
print("Method: Classical fourth-order Runge-Kutta")
print(f"Requested step size: {h}")
print(f"Actual step size: {h_actual:.6f}")

print("\nERROR ANALYSIS")
print("-" * 45)
print(f"Maximum error in x: {np.max(errors_x):.12e}")
print(f"Maximum error in y: {np.max(errors_y):.12e}")
print(f"Maximum energy error: {np.max(energy_error):.12e}")

print("\nSELECTED NUMERICAL RESULTS")
print("-" * 65)
print(f"{'t':>8} {'Numerical x':>16} {'Exact x':>16} "
      f"{'Numerical y':>16}")

for i in range(0, len(t), 100):
    print(
        f"{t[i]:8.2f} "
        f"{x_numerical[i]:16.10f} "
        f"{x_true[i]:16.10f} "
        f"{y_numerical[i]:16.10f}"
    )

# Include final point
i = len(t) - 1
print(
    f"{t[i]:8.2f} "
    f"{x_numerical[i]:16.10f} "
    f"{x_true[i]:16.10f} "
    f"{y_numerical[i]:16.10f}"
)


# Fourth-order convergence study
step_sizes = [0.2, 0.1, 0.05, 0.025]
max_errors = []

print("\nFOURTH-ORDER CONVERGENCE STUDY")
print("-" * 60)
print(f"{'Step size':>12} {'Maximum error':>20} {'Observed order':>18}")

for i, step in enumerate(step_sizes):

    t_test, z_test, _ = rk4(system, t0, z0, t_end, step)

    error = max(
        np.max(np.abs(z_test[:, 0] - exact_x(t_test))),
        np.max(np.abs(z_test[:, 1] - exact_y(t_test)))
    )

    max_errors.append(error)

    if i == 0:
        print(f"{step:12.5f} {error:20.8e} {'---':>18}")
    else:
        order = np.log(max_errors[i-1] / error) / np.log(2)
        print(f"{step:12.5f} {error:20.8e} {order:18.4f}")


# Plot 1: Numerical and exact solutions
plt.figure(figsize=(10, 5))

plt.plot(t, x_numerical, label="Numerical x(t)")
plt.plot(t, x_true, "--", label="Exact x(t)")
plt.plot(t, y_numerical, label="Numerical y(t)")
plt.plot(t, y_true, "--", label="Exact y(t)")

plt.xlabel("t")
plt.ylabel("Solutions")
plt.title("System of ODEs Using Classical RK4")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("systems_of_odes_rk4.png", dpi=300)
plt.show()


# Plot 2: Absolute errors
plt.figure(figsize=(9, 5))

plt.semilogy(t, errors_x, label="Error in x(t)")
plt.semilogy(t, errors_y, label="Error in y(t)")

plt.xlabel("t")
plt.ylabel("Absolute error (log scale)")
plt.title("RK4 Error Analysis")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("systems_of_odes_rk4_error.png", dpi=300)
plt.show()


# Plot 3: Phase-plane trajectory
plt.figure(figsize=(7, 7))

plt.plot(x_numerical, y_numerical, label="Numerical trajectory")
plt.plot(x_true, y_true, "--", label="Exact trajectory")

plt.scatter([z0[0]], [z0[1]], label="Initial point")

plt.xlabel("x(t)")
plt.ylabel("y(t)")
plt.title("Phase-Plane Trajectory Using RK4")
plt.axis("equal")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("systems_of_odes_rk4_phase_plane.png", dpi=300)
plt.show()


# Plot 4: Energy conservation
plt.figure(figsize=(9, 5))

plt.plot(t, energy, label="Numerical energy")
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
plt.savefig("systems_of_odes_rk4_energy.png", dpi=300)
plt.show()


# Plot 5: Fourth-order convergence
plt.figure(figsize=(9, 5))

plt.loglog(step_sizes, max_errors, "o-", label="RK4 maximum error")

reference = [
    max_errors[0] * (step / step_sizes[0])**4
    for step in step_sizes
]

plt.loglog(step_sizes, reference, "--", label="Fourth-order reference")

plt.xlabel("Step size")
plt.ylabel("Maximum absolute error")
plt.title("Fourth-Order Convergence of RK4")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()
plt.savefig("systems_of_odes_rk4_convergence.png", dpi=300)
plt.show()

