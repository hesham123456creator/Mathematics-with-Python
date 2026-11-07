
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# 1. Physical and numerical parameters
# --------------------------------------------------

L = 1.0
T = 0.1
alpha = 1.0

Nx = 50
dx = L / Nx

# Choose a stable time step
target_r = 0.4
dt_target = target_r * dx**2 / alpha

Nt = int(np.ceil(T / dt_target))
dt = T / Nt

r = alpha * dt / dx**2

if Nx < 2 or T <= 0 or alpha <= 0:
    raise ValueError("Invalid physical or numerical parameters.")

if r > 0.5:
    raise ValueError(
        f"Stability condition violated: r = {r:.4f}"
    )

print("HEAT DIFFUSION: EXPLICIT FINITE DIFFERENCE METHOD")
print("-" * 50)
print(f"Spatial steps: {Nx}")
print(f"Time steps: {Nt}")
print(f"dx = {dx:.6f}")
print(f"dt = {dt:.8f}")
print(f"Stability parameter r = {r:.6f}")

# --------------------------------------------------
# 2. Spatial and time grids
# --------------------------------------------------

x = np.linspace(0, L, Nx + 1)
time = np.linspace(0, T, Nt + 1)

# Initial condition
u = np.sin(np.pi * x)

# Boundary conditions
u[0] = 0.0
u[-1] = 0.0

# Selected time steps for visualization
selected_steps = sorted(set([
    0,
    Nt // 4,
    Nt // 2,
    3 * Nt // 4,
    Nt
]))

snapshots = {0: u.copy()}

# --------------------------------------------------
# 3. Exact solution
# --------------------------------------------------

def exact_solution(x, t):
    return (
        np.exp(-alpha * np.pi**2 * t)
        * np.sin(np.pi * x)
    )

# --------------------------------------------------
# 4. Numerical solution: FTCS method
# --------------------------------------------------

for n in range(Nt):

    u_new = u.copy()

    # Vectorized finite difference formula
    u_new[1:-1] = (
        u[1:-1]
        + r * (
            u[2:]
            - 2 * u[1:-1]
            + u[:-2]
        )
    )

    # Enforce boundary conditions
    u_new[0] = 0.0
    u_new[-1] = 0.0

    u = u_new

    if n + 1 in selected_steps:
        snapshots[n + 1] = u.copy()

# --------------------------------------------------
# 5. Error analysis
# --------------------------------------------------

u_exact = exact_solution(x, T)

absolute_error = np.abs(u - u_exact)

maximum_error = np.max(absolute_error)
mean_error = np.mean(absolute_error)
l2_error = np.sqrt(np.mean(absolute_error**2))

print("\nERROR ANALYSIS")
print("-" * 50)
print(f"Maximum absolute error: {maximum_error:.8e}")
print(f"Mean absolute error:    {mean_error:.8e}")
print(f"RMS error:              {l2_error:.8e}")

# --------------------------------------------------
# 6. Temperature evolution
# --------------------------------------------------

plt.figure(figsize=(10, 6))

for step in selected_steps:

    t = time[step]

    plt.plot(
        x,
        snapshots[step],
        label=f"Numerical, t={t:.3f}"
    )

    plt.plot(
        x,
        exact_solution(x, t),
        "--",
        label=f"Exact, t={t:.3f}"
    )

plt.xlabel("Position x")
plt.ylabel("Temperature u(x,t)")
plt.title("One-Dimensional Heat Diffusion")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("heat_diffusion_evolution.png", dpi=300)
plt.show()

# --------------------------------------------------
# 7. Final numerical error
# --------------------------------------------------

plt.figure(figsize=(9, 5))

plt.plot(x, absolute_error, marker="o", markersize=3)

plt.xlabel("Position x")
plt.ylabel("Absolute error")
plt.title("Numerical Error at Final Time")
plt.grid(True)
plt.tight_layout()
plt.savefig("heat_diffusion_error.png", dpi=300)
plt.show()

# --------------------------------------------------
# 8. Convergence study
# --------------------------------------------------

grid_sizes = [20, 40, 80, 160, 320]

errors = []
grid_spacings = []

print("\nCONVERGENCE STUDY")
print("-" * 65)
print(f"{'Nx':>8} {'dx':>12} {'dt':>12} {'Max Error':>16}")

for N in grid_sizes:

    h = L / N

    target_dt = target_r * h**2 / alpha
    steps = int(np.ceil(T / target_dt))
    step_size = T / steps

    ratio = alpha * step_size / h**2

    grid = np.linspace(0, L, N + 1)

    numerical = np.sin(np.pi * grid)

    for n in range(steps):

        updated = numerical.copy()

        updated[1:-1] = (
            numerical[1:-1]
            + ratio * (
                numerical[2:]
                - 2 * numerical[1:-1]
                + numerical[:-2]
            )
        )

        numerical = updated

    exact = exact_solution(grid, T)

    error = np.max(np.abs(numerical - exact))

    grid_spacings.append(h)
    errors.append(error)

    print(
        f"{N:8d} {h:12.6f} "
        f"{step_size:12.8f} {error:16.8e}"
    )

# Calculate observed convergence orders
print("\nOBSERVED CONVERGENCE ORDERS")

for i in range(1, len(errors)):

    order = (
        np.log(errors[i - 1] / errors[i])
        / np.log(grid_spacings[i - 1] / grid_spacings[i])
    )

    print(
        f"N = {grid_sizes[i - 1]:3d} -> "
        f"{grid_sizes[i]:3d}: order = {order:.4f}"
    )

# --------------------------------------------------
# 9. Convergence plot
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.loglog(
    grid_spacings,
    errors,
    "o-",
    label="Numerical error"
)

plt.xlabel("Grid spacing dx")
plt.ylabel("Maximum absolute error")
plt.title("Convergence of the FTCS Method")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig("heat_diffusion_convergence.png", dpi=300)
plt.show()

