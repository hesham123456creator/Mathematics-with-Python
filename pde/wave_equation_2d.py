import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --------------------------------------------------
# Physical parameters
# --------------------------------------------------

Lx = 1.0
Ly = 1.0
T = 0.2
c = 1.0

# Ensure output directory exists
Path("images").mkdir(exist_ok=True)

# --------------------------------------------------
# Numerical parameters
# --------------------------------------------------

Nx = 50
Ny = 50

dx = Lx / Nx
dy = Ly / Ny

# CFL stability condition:
# c^2 * dt^2 * (1/dx^2 + 1/dy^2) <= 1

safety_factor = 0.9

dt_limit = safety_factor / (
    c * np.sqrt(1 / dx**2 + 1 / dy**2)
)

Nt = int(np.ceil(T / dt_limit))
dt = T / Nt

rx2 = (c * dt / dx)**2
ry2 = (c * dt / dy)**2

stability_parameter = rx2 + ry2

print("Stability Analysis")
print("-" * 35)
print(f"rx^2 = {rx2:.10f}")
print(f"ry^2 = {ry2:.10f}")
print(f"Stability parameter = {stability_parameter:.10f}")
print("Stability condition: rx^2 + ry^2 <= 1")

# --------------------------------------------------
# Spatial grid
# --------------------------------------------------

x = np.linspace(0, Lx, Nx + 1)
y = np.linspace(0, Ly, Ny + 1)

X, Y = np.meshgrid(x, y)

# --------------------------------------------------
# Initial conditions
# --------------------------------------------------

u_initial = np.sin(np.pi * X) * np.sin(np.pi * Y)

velocity_initial = np.zeros_like(u_initial)

# --------------------------------------------------
# Initial time levels
# --------------------------------------------------

u_previous = u_initial.copy()
u_current = u_initial.copy()

# --------------------------------------------------
# First time step
# --------------------------------------------------

u_current[1:-1, 1:-1] = (
    u_initial[1:-1, 1:-1]
    + dt * velocity_initial[1:-1, 1:-1]
    + 0.5 * rx2 * (
        u_initial[1:-1, 2:]
        - 2 * u_initial[1:-1, 1:-1]
        + u_initial[1:-1, :-2]
    )
    + 0.5 * ry2 * (
        u_initial[2:, 1:-1]
        - 2 * u_initial[1:-1, 1:-1]
        + u_initial[:-2, 1:-1]
    )
)

# --------------------------------------------------
# Time evolution
# --------------------------------------------------

for n in range(1, Nt):

    u_next = np.zeros_like(u_current)

    laplacian_x = (
        u_current[1:-1, 2:]
        - 2 * u_current[1:-1, 1:-1]
        + u_current[1:-1, :-2]
    )

    laplacian_y = (
        u_current[2:, 1:-1]
        - 2 * u_current[1:-1, 1:-1]
        + u_current[:-2, 1:-1]
    )

    u_next[1:-1, 1:-1] = (
        2 * u_current[1:-1, 1:-1]
        - u_previous[1:-1, 1:-1]
        + rx2 * laplacian_x
        + ry2 * laplacian_y
    )

    u_previous = u_current.copy()
    u_current = u_next.copy()

# --------------------------------------------------
# Exact analytical solution
# --------------------------------------------------

u_exact = (
    np.sin(np.pi * X)
    * np.sin(np.pi * Y)
    * np.cos(np.sqrt(2) * np.pi * c * T)
)

# --------------------------------------------------
# Accuracy analysis
# --------------------------------------------------

error = u_current - u_exact
absolute_error = np.abs(error)

maximum_error = np.max(absolute_error)
mean_error = np.mean(absolute_error)
rms_error = np.sqrt(np.mean(error**2))

print("\nAccuracy Analysis")
print("-" * 35)
print(f"Maximum absolute error: {maximum_error:.10e}")
print(f"Mean absolute error: {mean_error:.10e}")
print(f"Root mean square error: {rms_error:.10e}")

# --------------------------------------------------
# Figure 1: Numerical solution
# --------------------------------------------------

fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection="3d")

surface = ax.plot_surface(
    X, Y, u_current,
    cmap="viridis",
    edgecolor="none"
)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("Displacement")
ax.set_title("2D Wave Equation: Numerical Solution")

fig.colorbar(surface, ax=ax, shrink=0.7)
plt.tight_layout()

plt.savefig("images/wave_equation_2d_solution.png", dpi=300)
plt.close()

# --------------------------------------------------
# Figure 2: Numerical versus exact solution
# --------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

contour1 = axes[0].contourf(
    X, Y, u_current, levels=30, cmap="viridis"
)

axes[0].set_title("Numerical Solution")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")

contour2 = axes[1].contourf(
    X, Y, u_exact, levels=30, cmap="viridis"
)

axes[1].set_title("Exact Solution")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")

fig.colorbar(contour1, ax=axes[0])
fig.colorbar(contour2, ax=axes[1])

plt.tight_layout()
plt.savefig("images/wave_equation_2d_comparison.png", dpi=300)
plt.close()

# --------------------------------------------------
# Figure 3: Absolute error
# --------------------------------------------------

plt.figure(figsize=(8, 6))

contour = plt.contourf(
    X, Y, absolute_error,
    levels=30,
    cmap="inferno"
)

plt.colorbar(contour, label="Absolute error")
plt.xlabel("x")
plt.ylabel("y")
plt.title("2D Wave Equation: Absolute Error")
plt.tight_layout()

plt.savefig("images/wave_equation_2d_error.png", dpi=300)
plt.close()

print("\nGenerated figures:")
print("1. images/wave_equation_2d_solution.png")
print("2. images/wave_equation_2d_comparison.png")
print("3. images/wave_equation_2d_error.png")

# --------------------------------------------------
# Convergence study
# --------------------------------------------------

def solve_wave_2d(N):

    dx_local = Lx / N
    dy_local = Ly / N

    x_local = np.linspace(0, Lx, N + 1)
    y_local = np.linspace(0, Ly, N + 1)

    X_local, Y_local = np.meshgrid(x_local, y_local)

    dt_max = safety_factor / (
        c * np.sqrt(1 / dx_local**2 + 1 / dy_local**2)
    )

    Nt_local = int(np.ceil(T / dt_max))
    dt_local = T / Nt_local

    rx_local = (c * dt_local / dx_local)**2
    ry_local = (c * dt_local / dy_local)**2

    u0 = (
        np.sin(np.pi * X_local)
        * np.sin(np.pi * Y_local)
    )

    velocity = np.zeros_like(u0)

    u_previous_local = u0.copy()
    u_current_local = u0.copy()

    # First time step
    u_current_local[1:-1, 1:-1] = (
        u0[1:-1, 1:-1]
        + dt_local * velocity[1:-1, 1:-1]
        + 0.5 * rx_local * (
            u0[1:-1, 2:]
            - 2 * u0[1:-1, 1:-1]
            + u0[1:-1, :-2]
        )
        + 0.5 * ry_local * (
            u0[2:, 1:-1]
            - 2 * u0[1:-1, 1:-1]
            + u0[:-2, 1:-1]
        )
    )

    # Time evolution
    for n in range(1, Nt_local):

        u_next_local = np.zeros_like(u_current_local)

        lap_x = (
            u_current_local[1:-1, 2:]
            - 2 * u_current_local[1:-1, 1:-1]
            + u_current_local[1:-1, :-2]
        )

        lap_y = (
            u_current_local[2:, 1:-1]
            - 2 * u_current_local[1:-1, 1:-1]
            + u_current_local[:-2, 1:-1]
        )

        u_next_local[1:-1, 1:-1] = (
            2 * u_current_local[1:-1, 1:-1]
            - u_previous_local[1:-1, 1:-1]
            + rx_local * lap_x
            + ry_local * lap_y
        )

        u_previous_local = u_current_local.copy()
        u_current_local = u_next_local.copy()

    exact_local = (
        np.sin(np.pi * X_local)
        * np.sin(np.pi * Y_local)
        * np.cos(np.sqrt(2) * np.pi * c * T)
    )

    return np.max(np.abs(u_current_local - exact_local))


# --------------------------------------------------
# Grid refinement results
# --------------------------------------------------

grid_sizes = [20, 40, 80, 160]

errors = []
spacings = []

print("\nConvergence Study")
print("=" * 65)

for N in grid_sizes:

    dx_local = Lx / N
    error_local = solve_wave_2d(N)

    spacings.append(dx_local)
    errors.append(error_local)

print("\nGrid Refinement Results")
print("-" * 65)
print(f"{'N':>8} {'dx':>12} {'Maximum Error':>18} {'Order':>12}")

for j in range(len(grid_sizes)):

    if j == 0:
        order = "---"
    else:
        order_value = (
            np.log(errors[j - 1] / errors[j])
            / np.log(spacings[j - 1] / spacings[j])
        )
        order = f"{order_value:.4f}"

    print(
        f"{grid_sizes[j]:8d} "
        f"{spacings[j]:12.6f} "
        f"{errors[j]:18.10e} "
        f"{order:>12}"
    )

# --------------------------------------------------
# Figure 4: Convergence plot
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.loglog(
    spacings,
    errors,
    "o-",
    linewidth=2,
    markersize=7,
    label="Numerical error"
)

reference_error = errors[-1] * (
    np.array(spacings) / spacings[-1]
)**2

plt.loglog(
    spacings,
    reference_error,
    "--",
    linewidth=2,
    label="Second-order reference"
)

plt.xlabel("Spatial step size")
plt.ylabel("Maximum absolute error")
plt.title("Convergence Study: 2D Wave Equation")
plt.legend()
plt.grid(True, which="both")
plt.tight_layout()

plt.savefig("images/wave_equation_2d_convergence.png", dpi=300)
plt.close()

print("\nGenerated figure:")
print("images/wave_equation_2d_convergence.png")
