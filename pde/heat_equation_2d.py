
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# Two-Dimensional Heat Equation
#
# u_t = alpha * (u_xx + u_yy)
#
# Initial condition:
# u(x, y, 0) = sin(pi*x) sin(pi*y)
#
# Boundary conditions:
# u = 0 on all four boundaries
#
# Exact solution:
# u(x, y, t) = exp(-2*alpha*pi^2*t)
#              * sin(pi*x) * sin(pi*y)
# --------------------------------------------------

# Physical parameters
Lx = 1.0
Ly = 1.0
alpha = 1.0
T = 0.05

# Numerical parameters
Nx = 50
Ny = 50

dx = Lx / Nx
dy = Ly / Ny

# Choose time step to satisfy stability condition
dt = 0.9 / (
    2.0 * alpha * (1.0 / dx**2 + 1.0 / dy**2)
)

Nt = int(np.ceil(T / dt))
dt = T / Nt

rx = alpha * dt / dx**2
ry = alpha * dt / dy**2

print("Stability Analysis")
print("------------------------------")
print(f"rx = {rx:.10f}")
print(f"ry = {ry:.10f}")
print(f"rx + ry = {rx + ry:.10f}")

if rx + ry > 0.5:
    raise ValueError("Stability condition violated.")

# Spatial grid
x = np.linspace(0, Lx, Nx + 1)
y = np.linspace(0, Ly, Ny + 1)

X, Y = np.meshgrid(x, y, indexing="ij")

# Initial condition
u = np.sin(np.pi * X) * np.sin(np.pi * Y)

# Exact analytical solution
def exact_solution(X, Y, t):
    return (
        np.exp(-2 * alpha * np.pi**2 * t)
        * np.sin(np.pi * X)
        * np.sin(np.pi * Y)
    )

# Time integration
for n in range(Nt):
    u_new = u.copy()

    u_new[1:-1, 1:-1] = (
        u[1:-1, 1:-1]
        + rx * (
            u[2:, 1:-1]
            - 2 * u[1:-1, 1:-1]
            + u[:-2, 1:-1]
        )
        + ry * (
            u[1:-1, 2:]
            - 2 * u[1:-1, 1:-1]
            + u[1:-1, :-2]
        )
    )

    u = u_new

# Exact solution at final time
u_exact = exact_solution(X, Y, T)

# Accuracy analysis
error = np.abs(u - u_exact)

maximum_error = np.max(error)
mean_error = np.mean(error)
rms_error = np.sqrt(np.mean(error**2))

print("\nAccuracy Analysis")
print("------------------------------")
print(f"Maximum absolute error: {maximum_error:.10e}")
print(f"Mean absolute error: {mean_error:.10e}")
print(f"Root mean square error: {rms_error:.10e}")

# --------------------------------------------------
# Visualization 1: Numerical solution
# --------------------------------------------------

plt.figure(figsize=(9, 7))

plt.contourf(X, Y, u, levels=30, cmap="hot")
plt.colorbar(label="Temperature")

plt.xlabel("x")
plt.ylabel("y")
plt.title(f"2D Heat Equation: Numerical Solution at t={T}")

plt.tight_layout()
plt.savefig("images/heat_equation_2d_solution.png", dpi=300)
plt.close()

# --------------------------------------------------
# Visualization 2: Numerical vs exact solution
# --------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

im1 = axes[0].contourf(X, Y, u, levels=30, cmap="hot")
axes[0].set_title("Numerical Solution")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")
fig.colorbar(im1, ax=axes[0])

im2 = axes[1].contourf(X, Y, u_exact, levels=30, cmap="hot")
axes[1].set_title("Exact Solution")
axes[1].set_xlabel("x")
axes[1].set_ylabel("y")
fig.colorbar(im2, ax=axes[1])

plt.tight_layout()
plt.savefig("images/heat_equation_2d_comparison.png", dpi=300)
plt.close()

# --------------------------------------------------
# Visualization 3: Absolute error
# --------------------------------------------------

plt.figure(figsize=(9, 7))

plt.contourf(X, Y, error, levels=30, cmap="viridis")
plt.colorbar(label="Absolute Error")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Absolute Numerical Error")

plt.tight_layout()
plt.savefig("images/heat_equation_2d_error.png", dpi=300)
plt.close()

print("\nGenerated figures:")
print("1. images/heat_equation_2d_solution.png")
print("2. images/heat_equation_2d_comparison.png")
print("3. images/heat_equation_2d_error.png")

# --------------------------------------------------
# Grid Refinement and Convergence Study
# --------------------------------------------------

def solve_heat_2d(N):
    dx = Lx / N
    dy = Ly / N

    dt_max = 0.9 / (
        2 * alpha * (1 / dx**2 + 1 / dy**2)
    )

    Nt = int(np.ceil(T / dt_max))
    dt = T / Nt

    rx = alpha * dt / dx**2
    ry = alpha * dt / dy**2

    x = np.linspace(0, Lx, N + 1)
    y = np.linspace(0, Ly, N + 1)

    X, Y = np.meshgrid(x, y, indexing="ij")

    u = np.sin(np.pi * X) * np.sin(np.pi * Y)

    for n in range(Nt):
        u_new = u.copy()

        u_new[1:-1, 1:-1] = (
            u[1:-1, 1:-1]
            + rx * (
                u[2:, 1:-1]
                - 2 * u[1:-1, 1:-1]
                + u[:-2, 1:-1]
            )
            + ry * (
                u[1:-1, 2:]
                - 2 * u[1:-1, 1:-1]
                + u[1:-1, :-2]
            )
        )

        u = u_new

    u_exact = exact_solution(X, Y, T)

    maximum_error = np.max(np.abs(u - u_exact))

    return dx, maximum_error


print("\nConvergence Study")
print("=" * 65)
print(f"{'N':>8} {'dx':>12} {'Maximum Error':>20} {'Order':>12}")
print("-" * 65)

grid_sizes = [25, 50, 100, 200]

errors = []
grid_spacings = []

for N in grid_sizes:
    dx, error = solve_heat_2d(N)

    grid_spacings.append(dx)
    errors.append(error)

    if len(errors) == 1:
        order_text = "---"
    else:
        order = np.log(errors[-2] / errors[-1]) / np.log(2)
        order_text = f"{order:.4f}"

    print(
        f"{N:8d} {dx:12.6f} "
        f"{error:20.10e} {order_text:>12}"
    )

# --------------------------------------------------
# Convergence Plot
# --------------------------------------------------

plt.figure(figsize=(8, 6))

plt.loglog(
    grid_spacings,
    errors,
    "o-",
    label="Numerical error"
)

reference = (
    errors[0]
    * (np.array(grid_spacings) / grid_spacings[0])**2
)

plt.loglog(
    grid_spacings,
    reference,
    "--",
    label="Second-order reference"
)

plt.xlabel("Grid spacing (dx)")
plt.ylabel("Maximum absolute error")
plt.title("2D Heat Equation: Grid Convergence")
plt.grid(True, which="both")
plt.legend()

plt.tight_layout()
plt.savefig(
    "images/heat_equation_2d_convergence.png",
    dpi=300
)
plt.close()

print("\nGenerated figure:")
print("images/heat_equation_2d_convergence.png")
