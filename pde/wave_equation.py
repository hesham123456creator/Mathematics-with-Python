
import numpy as np
import matplotlib.pyplot as plt

# Physical parameters
L = 1.0
T = 1.0
c = 1.0

# Numerical parameters
Nx = 50
dx = L / Nx

# Choose time step to satisfy the stability condition
dt = 0.9 * dx / c
Nt = int(T / dt)
dt = T / Nt

r = c * dt / dx

print("Stability parameter r =", r)

# Spatial grid
x = np.linspace(0, L, Nx + 1)

# Initial displacement and velocity
u_initial = np.sin(np.pi * x)
velocity_initial = np.zeros(Nx + 1)

# Initial time levels
u_previous = u_initial.copy()
u_current = u_initial.copy()

# Apply the first time step
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

# Time evolution
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

# Exact solution
u_exact = (
    np.cos(np.pi * c * T)
    * np.sin(np.pi * x)
)

# Calculate maximum error
error = np.max(np.abs(u_current - u_exact))

print("Maximum absolute error =", error)

# Plot numerical and exact solutions
plt.figure(figsize=(9, 6))

plt.plot(x, u_current, "o", label="Numerical solution")
plt.plot(x, u_exact, label="Exact solution")

plt.xlabel("Position x")
plt.ylabel("Displacement u(x,T)")
plt.title("One-Dimensional Wave Equation")
plt.grid(True)
plt.legend()
plt.show()
