# Pylint does not fully recognize SciPy's dynamically exposed functions.
# pylint: disable=no-name-in-module
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad, trapezoid, simpson
from scipy.special import erf

# --------------------------------------------------
# 1. Define the function and integration interval
# --------------------------------------------------

def f(x):
    return np.exp(-x**2)


a = 0.0
b = 1.0

# Exact value using the error function
exact_value = np.sqrt(np.pi) / 2 * erf(1)

print("=" * 65)
print("NUMERICAL INTEGRATION: ACCURACY AND CONVERGENCE")
print("=" * 65)

print(f"\nIntegration interval: [{a}, {b}]")
print(f"Exact integral: {exact_value:.12f}")

# --------------------------------------------------
# 2. Adaptive quadrature using SciPy
# --------------------------------------------------

quad_result, quad_error = quad(
    f,
    a,
    b,
    epsabs=1e-12,
    epsrel=1e-12
)

print("\nADAPTIVE QUADRATURE")
print("-" * 65)
print(f"Numerical result: {quad_result:.12f}")
print(f"Estimated error:  {quad_error:.8e}")
print(f"Actual error:     {abs(quad_result - exact_value):.8e}")

# --------------------------------------------------
# 3. Composite trapezoidal and Simpson's rules
# --------------------------------------------------

N = 100

x = np.linspace(a, b, N + 1)
y = f(x)

trapezoidal_result = trapezoid(y, x=x)
simpson_result = simpson(y, x=x)

print("\nCOMPOSITE NUMERICAL METHODS")
print("-" * 65)

print(f"Number of subintervals: {N}")

print(f"\nTrapezoidal rule: {trapezoidal_result:.12f}")
print(
    "Absolute error:",
    f"{abs(trapezoidal_result - exact_value):.8e}"
)

print(f"\nSimpson's rule:   {simpson_result:.12f}")
print(
    "Absolute error:",
    f"{abs(simpson_result - exact_value):.8e}"
)

# --------------------------------------------------
# 4. Convergence study
# --------------------------------------------------

interval_counts = [10, 20, 40, 80, 160, 320, 640]

trapezoidal_errors = []
simpson_errors = []

print("\nCONVERGENCE STUDY")
print("-" * 80)

print(
    f"{'N':>8} {'Trapezoidal Error':>22} "
    f"{'Simpson Error':>22}"
)

for n in interval_counts:

    grid = np.linspace(a, b, n + 1)
    values = f(grid)

    trap_result = trapezoid(values, x=grid)
    simp_result = simpson(values, x=grid)

    trap_error = abs(trap_result - exact_value)
    simp_error = abs(simp_result - exact_value)

    trapezoidal_errors.append(trap_error)
    simpson_errors.append(simp_error)

    print(
        f"{n:8d} {trap_error:22.8e} "
        f"{simp_error:22.8e}"
    )

# --------------------------------------------------
# 5. Observed convergence orders
# --------------------------------------------------

print("\nOBSERVED CONVERGENCE ORDERS")
print("-" * 65)

for i in range(1, len(interval_counts)):

    trap_order = np.log(
        trapezoidal_errors[i - 1]
        / trapezoidal_errors[i]
    ) / np.log(2)

    simp_order = np.log(
        simpson_errors[i - 1]
        / simpson_errors[i]
    ) / np.log(2)

    print(
        f"N = {interval_counts[i - 1]:3d} -> "
        f"{interval_counts[i]:3d}: "
        f"Trapezoidal = {trap_order:.4f}, "
        f"Simpson = {simp_order:.4f}"
    )

# --------------------------------------------------
# 6. Plot the integrand
# --------------------------------------------------

x_plot = np.linspace(a, b, 400)

plt.figure(figsize=(9, 5))

plt.plot(
    x_plot,
    f(x_plot),
    linewidth=2,
    label=r"$f(x)=e^{-x^2}$"
)

plt.fill_between(
    x_plot,
    f(x_plot),
    alpha=0.2,
    label="Area under the curve"
)

plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Numerical Integration of the Gaussian Function")
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig(
    "numerical_integration_function.png",
    dpi=300
)

plt.show()

# --------------------------------------------------
# 7. Plot convergence
# --------------------------------------------------

plt.figure(figsize=(9, 5))

plt.loglog(
    interval_counts,
    trapezoidal_errors,
    "o-",
    label="Trapezoidal rule"
)

plt.loglog(
    interval_counts,
    simpson_errors,
    "s-",
    label="Simpson's rule"
)

plt.xlabel("Number of subintervals")
plt.ylabel("Absolute error")
plt.title("Convergence of Numerical Integration Methods")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()

plt.savefig(
    "numerical_integration_convergence.png",
    dpi=300
)

plt.show()

print("\nGenerated figures:")
print("- numerical_integration_function.png")
print("- numerical_integration_convergence.png")

print("\nNumerical integration completed successfully.")
