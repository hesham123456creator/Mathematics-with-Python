import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# Physical parameters
# --------------------------------------------------

L = 1.0       # Length of the spatial domain
T = 1.0       # Final simulation time
c = 1.0       # Wave propagation speed

# --------------------------------------------------
# Numerical parameters
# --------------------------------------------------

Nx = 50                  # Number of spatial intervals
dx = L / Nx              # Spatial step size

# Choose time step to satisfy the CFL stability condition
dt = 0.9 * dx / c

Nt = int(T / dt)          # Number of time steps
dt = T / Nt               # Adjust time step to reach exactly T

# CFL stability parameter
r = c * dt / dx

print("Stability parameter r =", r)

# --------------------------------------------------
# Spatial grid
# --------------------------------------------------

x = np.linspace(0, L, Nx + 1)

# --------------------------------------------------
# Initial conditions
# --------------------------------------------------

# Initial displacement: u(x, 0) = sin(pi*x)
u_initial = np.sin(np.pi * x)

# Initial velocity: u_t(x, 0) = 0
velocity_initial = np.zeros(Nx + 1)
# --------------------------------------------------
# Initial time levels
# --------------------------------------------------

u_previous = u_initial.copy()
u_current = u_initial.copy()

# --------------------------------------------------
# First time step
# --------------------------------------------------

for i in range(1, Nx):
    u_current[i] = (
        u_initial[i]
        + dt * velocity_initial[i]
        + 0.5 * r**2 * (
            u_initial[i + 1]
            - 2 * u_initial[i]
            + u_initial[i - 1]
        )
    )

# --------------------------------------------------
# Time evolution
# --------------------------------------------------

for n in range(1, Nt):

    u_next = np.zeros(Nx + 1)

    for i in range(1, Nx):
        u_next[i] = (
            2 * u_current[i]
            - u_previous[i]
            + r**2 * (
                u_current[i + 1]
                - 2 * u_current[i]
                + u_current[i - 1]
            )
        )

    u_previous = u_current.copy()
    u_current = u_next.copy()

# --------------------------------------------------
# Exact analytical solution
# --------------------------------------------------

u_exact = np.sin(np.pi * x) * np.cos(np.pi * c * T)

# --------------------------------------------------
# Accuracy analysis
# --------------------------------------------------

error = u_current - u_exact

maximum_error = np.max(np.abs(error))
rms_error = np.sqrt(np.mean(error**2))

print("\nAccuracy Analysis")
print("------------------------------")
print(f"Maximum absolute error: {maximum_error:.10e}")
print(f"Root mean square error: {rms_error:.10e}")

# --------------------------------------------------
# Plot numerical and exact solutions
# --------------------------------------------------

plt.figure(figsize=(9, 5))

plt.plot(
    x, u_current,
    "o",
    markersize=4,
    label="Numerical solution"
)

plt.plot(
    x, u_exact,
    "-",
    linewidth=2,
    label="Exact solution"
)

plt.xlabel("Position x")
plt.ylabel("Displacement u(x,T)")
plt.title("One-Dimensional Wave Equation")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("images/wave_equation_comparison.png", dpi=300)
plt.show()

# --------------------------------------------------
# Convergence study
# --------------------------------------------------

T_conv = 0.5

Nx_values = [25, 50, 100, 200, 400]

errors = []
grid_sizes = []

print("\nConvergence Study")
print("=" * 60)
for Nx_conv in Nx_values:

    dx_conv = L / Nx_conv

    dt_conv = 0.9 * dx_conv / c
    Nt_conv = int(T_conv / dt_conv)
    dt_conv = T_conv / Nt_conv

    r_conv = c * dt_conv / dx_conv

    x_conv = np.linspace(0, L, Nx_conv + 1)

    # Initial conditions
    u_initial_conv = np.sin(np.pi * x_conv)
    velocity_conv = np.zeros(Nx_conv + 1)

    # Initial time levels
    u_previous_conv = u_initial_conv.copy()
    u_current_conv = u_initial_conv.copy()

    # First time step
    for i in range(1, Nx_conv):
        u_current_conv[i] = (
            u_initial_conv[i]
            + dt_conv * velocity_conv[i]
            + 0.5 * r_conv**2 * (
                u_initial_conv[i + 1]
                - 2 * u_initial_conv[i]
                + u_initial_conv[i - 1]
            )
        )

    # Time evolution
    for n in range(1, Nt_conv):

        u_next_conv = np.zeros(Nx_conv + 1)

        for i in range(1, Nx_conv):
            u_next_conv[i] = (
                2 * u_current_conv[i]
                - u_previous_conv[i]
                + r_conv**2 * (
                    u_current_conv[i + 1]
                    - 2 * u_current_conv[i]
                    + u_current_conv[i - 1]
                )
            )

        u_previous_conv = u_current_conv.copy()
        u_current_conv = u_next_conv.copy()

    # Exact solution
    u_exact_conv = (
        np.sin(np.pi * x_conv)
        * np.cos(np.pi * c * T_conv)
    )

    # Maximum absolute error
    error_conv = np.max(
        np.abs(u_current_conv - u_exact_conv)
    )

    errors.append(error_conv)
    grid_sizes.append(dx_conv)

# --------------------------------------------------
# Convergence order
# --------------------------------------------------

print("\nGrid Refinement Results")
print("-" * 65)
print(
    f"{'Nx':>8} {'dx':>12} "
    f"{'Maximum Error':>18} {'Order':>12}"
)

for j in range(len(Nx_values)):

    if j == 0:
        order = "---"
    else:
        order_value = (
            np.log(errors[j - 1] / errors[j])
            / np.log(grid_sizes[j - 1] / grid_sizes[j])
        )
        order = f"{order_value:.4f}"

    print(
        f"{Nx_values[j]:8d} "
        f"{grid_sizes[j]:12.6f} "
        f"{errors[j]:18.10e} "
        f"{order:>12}"
    )

# --------------------------------------------------
# Convergence plot
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.loglog(
    grid_sizes,
    errors,
    "o-",
    linewidth=2,
    markersize=7,
    label="Numerical error"
)

# Reference second-order slope
reference_error = errors[-1] * (
    np.array(grid_sizes) / grid_sizes[-1]
)**2

plt.loglog(
    grid_sizes,
    reference_error,
    "--",
    linewidth=2,
    label="Second-order reference"
)

plt.xlabel("Spatial step size (dx)")
plt.ylabel("Maximum absolute error")
plt.title("Convergence Study: Wave Equation")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()

plt.savefig(
    "images/wave_equation_convergence.png",
    dpi=300
)

plt.show()
