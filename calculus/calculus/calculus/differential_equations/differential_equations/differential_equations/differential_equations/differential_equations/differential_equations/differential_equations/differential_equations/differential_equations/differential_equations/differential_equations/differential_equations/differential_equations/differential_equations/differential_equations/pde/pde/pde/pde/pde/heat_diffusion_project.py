
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------
# 1. Physical parameters
# --------------------------------

L = 1.0
T = 0.1
alpha = 1.0

# --------------------------------
# 2. Numerical parameters
# --------------------------------

Nx = 50
dx = L / Nx

# Choose a stable time step
r = 0.4
dt = r * dx**2 / alpha

Nt = int(round(T / dt))
dt = T / Nt

r = alpha * dt / dx**2

print("Time step:", dt)
print("Number of time steps:", Nt)
print("Stability parameter:", r)

# --------------------------------
# 3. Spatial and time grids
# --------------------------------

x = np.linspace(0, L, Nx + 1)
time = np.linspace(0, T, Nt + 1)

# Initial temperature
u = np.sin(np.pi * x)

# Boundary conditions
u[0] = 0.0
u[-1] = 0.0

# Store solutions at selected times
selected_steps = {
    0,
    Nt // 4,
    Nt // 2,
    Nt
}

snapshots = {0: u.copy()}

# --------------------------------
# 4. Numerical solution
# --------------------------------

for n in range(Nt):

    u_new = u.copy()

    for i in range(1, Nx):

        u_new[i] = (
            u[i]
            + r * (
                u[i + 1]
                - 2 * u[i]
                + u[i - 1]
            )
        )

    u = u_new

    if n + 1 in selected_steps:
        snapshots[n + 1] = u.copy()

# --------------------------------
# 5. Exact solution
# --------------------------------

def exact_solution(x, t):
    return (
        np.exp(-alpha * np.pi**2 * t)
        * np.sin(np.pi * x)
    )

u_exact = exact_solution(x, T)

# --------------------------------
# 6. Error analysis
# --------------------------------

absolute_error = np.abs(u - u_exact)
maximum_error = np.max(absolute_error)
mean_error = np.mean(absolute_error)

print("\nError Analysis")
print("Maximum absolute error:", maximum_error)
print("Mean absolute error:", mean_error)

# --------------------------------
# 7. Plot temperature evolution
# --------------------------------

plt.figure(figsize=(10, 6))

for step in sorted(snapshots):

    t = step * dt

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
plt.ylabel("Temperature")
plt.title("Heat Diffusion Along a Rod")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# --------------------------------
# 8. Plot final numerical error
# --------------------------------

plt.figure(figsize=(9, 5))

plt.plot(x, absolute_error, marker="o")

plt.xlabel("Position x")
plt.ylabel("Absolute error")
plt.title("Numerical Error at Final Time")
plt.grid(True)
plt.tight_layout()
plt.show()
