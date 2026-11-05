
import numpy as np
import matplotlib.pyplot as plt

# Domain parameters
Lx = 1.0
Ly = 1.0

Nx = 30
Ny = 30

dx = Lx / Nx
dy = Ly / Ny

# Create the grid
x = np.linspace(0, Lx, Nx + 1)
y = np.linspace(0, Ly, Ny + 1)

X, Y = np.meshgrid(x, y)

# Initialize temperature
u = np.zeros((Ny + 1, Nx + 1))

# Boundary conditions
u[:, 0] = 0
u[:, -1] = 0
u[0, :] = 0
u[-1, :] = np.sin(np.pi * x)

# Iterative parameters
tolerance = 1e-6
max_iterations = 20000

# Gauss-Seidel iteration
for iteration in range(max_iterations):

    max_change = 0.0

    for j in range(1, Ny):
        for i in range(1, Nx):

            old_value = u[j, i]

            u[j, i] = (
                dy**2 * (u[j, i+1] + u[j, i-1])
                + dx**2 * (u[j+1, i] + u[j-1, i])
            ) / (2 * (dx**2 + dy**2))

            change = abs(u[j, i] - old_value)
            max_change = max(max_change, change)

    if max_change < tolerance:
        break

print("Number of iterations:", iteration + 1)
print("Maximum change:", max_change)

# Plot the numerical solution
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection="3d")

ax.plot_surface(X, Y, u, cmap="viridis")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("Temperature")
ax.set_title("Two-Dimensional Laplace Equation")

plt.show()

# Contour plot
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, u, levels=20, cmap="viridis")
plt.colorbar(label="Temperature")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Temperature Distribution")
plt.show()
