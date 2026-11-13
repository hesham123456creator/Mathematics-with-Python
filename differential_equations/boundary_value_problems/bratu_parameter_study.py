import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.optimize import brentq

# --------------------------------------------------
# Bratu equation: parameter study
#
# y'' + lambda * exp(y) = 0
# y(0) = y(1) = 0
# --------------------------------------------------

# Critical parameter
critical_a = brentq(
    lambda a: a * np.tanh(a / 2) - 2,
    1.0,
    5.0
)

critical_lambda = (
    2 * critical_a**2 / np.cosh(critical_a / 2)**2
)

print(f"Critical lambda: {critical_lambda:.10f}")

# Values below the critical parameter
lambda_values = [
    1.0,
    2.0,
    3.0,
    3.3,
    3.4,
    3.45,
    3.48,
    3.50,
    3.51,
    3.513,
    3.5135
]

x = np.linspace(0, 1, 100)
x_plot = np.linspace(0, 1, 300)

results = []

plt.figure(figsize=(10, 7))


def analytical_parameter(a):
    return 2 * a**2 / np.cosh(a / 2)**2


def exact_solution(x, a):
    return 2 * np.log(
        np.cosh(a / 2) /
        np.cosh(a * (x - 0.5))
    )


def exact_derivative(x, a):
    return -2 * a * np.tanh(a * (x - 0.5))


# --------------------------------------------------
# Study both solution branches
# --------------------------------------------------

for lam in lambda_values:

    # Find the two analytical parameters
    equation = lambda a: analytical_parameter(a) - lam

    a_lower = brentq(equation, 1e-8, critical_a)
    a_upper = brentq(equation, critical_a, 50)

    for branch, a in [
        ("Lower", a_lower),
        ("Upper", a_upper)
    ]:

        
        # Exact solution and initial guess
        y_exact_mesh = exact_solution(x, a)
        yp_exact_mesh = exact_derivative(x, a)

        guess_factor = 0.8 if branch == "Lower" else 1.2

        y_guess = np.vstack((
            guess_factor * y_exact_mesh,
            guess_factor * yp_exact_mesh
        ))

        # First-order system
        def ode(x, y):
            return np.vstack((
                y[1],
                -lam * np.exp(y[0])
            ))

        # Boundary conditions
        def bc(ya, yb):
            return np.array([
                ya[0],
                yb[0]
            ])

        # Numerical solution
        solution = solve_bvp(
            ode,
            bc,
            x,
            y_guess,
            tol=1e-8,
            max_nodes=2000
        )

        if not solution.success:
            raise RuntimeError(
                f"Solver failed for lambda={lam}, "
                f"branch={branch}: {solution.message}"
            )

        # Evaluate solutions
        y_numerical = solution.sol(x_plot)[0]
        y_exact = exact_solution(x_plot, a)

        maximum_error = np.max(
            np.abs(y_numerical - y_exact)
        )

        # ODE residual
        derivative = solution.sol(x_plot, 1)

        ode_residual = np.max(
            np.abs(
                derivative[1] + lam * np.exp(y_numerical)
            )
        )

        scipy_rms = np.max(solution.rms_residuals)

        results.append((
            lam,
            branch,
            a,
            np.max(y_numerical),
            maximum_error,
            ode_residual,
            scipy_rms,
            solution.x.size
        ))

        # Plot
        line_style = "-" if branch == "Lower" else "--"

        plt.plot(
            x_plot,
            y_numerical,
            linestyle=line_style,
            linewidth=2,
            label=f"lambda={lam}, {branch}"
        )

# --------------------------------------------------
# Display results
# --------------------------------------------------

print("\nBratu Parameter Study")
print("=" * 105)

print(
    f"{'Lambda':>10} {'Branch':>9} {'a':>10} "
    f"{'Maximum y':>12} {'Max Error':>14} "
    f"{'ODE Residual':>14} {'RMS Residual':>14} {'Nodes':>7}"
)

for row in results:
    lam, branch, a, ymax, error, residual, rms, nodes = row

    print(
        f"{lam:10.5f} {branch:>9} {a:10.6f} "
        f"{ymax:12.6f} {error:14.3e} "
        f"{residual:14.3e} {rms:14.3e} {nodes:7d}"
    )

# --------------------------------------------------
# Plot
# --------------------------------------------------

plt.xlabel("x")
plt.ylabel("y(x)")
plt.title("Bratu Problem: Lower and Upper Solution Branches")
plt.legend(ncol=2)
plt.grid(True)
plt.tight_layout()

plt.savefig(
    "images/bratu_bifurcation_study.png",
    dpi=300
)

plt.show()
plt.close()
# --------------------------------------------------
# Bifurcation Diagram
# --------------------------------------------------

plt.figure(figsize=(9, 6))

for branch in ["Lower", "Upper"]:

    branch_results = [
        row for row in results if row[1] == branch
    ]

    lambda_values = [row[0] for row in branch_results]
    maximum_values = [row[3] for row in branch_results]

    plt.plot(
        lambda_values,
        maximum_values,
        marker="o",
        linewidth=2,
        label=f"{branch} branch"
    )

# Critical point
critical_y = 2 * np.log(np.cosh(critical_a / 2))

plt.plot(
    critical_lambda,
    critical_y,
    "r*",
    markersize=14,
    label="Critical point"
)

plt.axvline(
    critical_lambda,
    color="gray",
    linestyle="--",
    alpha=0.7
)

plt.xlabel("Parameter lambda")
plt.ylabel("Maximum solution value")
plt.title("Bratu Problem: Bifurcation Diagram")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
    "images/bratu_bifurcation_diagram.png",
    dpi=300
)

plt.show()
plt.close()
