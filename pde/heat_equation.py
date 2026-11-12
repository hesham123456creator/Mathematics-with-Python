
import numpy as np
import matplotlib.pyplot as plt

# Physical and numerical parameters
L = 1.0
T = 0.1
alpha = 1.0

Nx = 50
Nt = 1000

dx = L / Nx
dt = T / Nt

# Stability parameter
r = alpha * dt / dx**2

print("Stability parameter r =", r)

if r > 0.5:
    raise ValueError("The numerical scheme is unstable. Reduce dt.")

# Spatial grid
x = np.linspace(0, L, Nx + 1)

# Initial temperature distribution
u = np.sin(np.pi * x)

# Boundary conditions
u[0] = 0
u[-1] = 0

# Time evolution using explicit finite differences
for n in range(Nt):
    u_new = u.copy()

    for i in range(1, Nx):
        u_new[i] = (
            u[i]
            + r * (u[i + 1] - 2 * u[i] + u[i - 1])
        )

    u = u_new

# Exact solution at final time
u_exact = (
    np.exp(-alpha * np.pi**2 * T)
    * np.sin(np.pi * x)
)

# Calculate maximum error
error = np.max(np.abs(u - u_exact))

print("Maximum absolute error =", error)

# Plot numerical and exact solutions
plt.figure(figsize=(9, 6))

plt.plot(x, u, "o", label="Numerical solution")
plt.plot(x, u_exact, label="Exact solution")

plt.xlabel("Position x")
plt.ylabel("Temperature u(x, T)")
plt.title("One-Dimensional Heat Equation")
plt.grid(True)
plt.legend()
plt.show()
