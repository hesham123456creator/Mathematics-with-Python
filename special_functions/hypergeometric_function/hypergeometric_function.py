import numpy as np
import matplotlib.pyplot as plt
from scipy.special import hyp2f1
from scipy.integrate import quad
import math


# ==================================================
# Gauss Hypergeometric Function
#
# 2F1(a,b;c;z)
# ==================================================


# --------------------------------------------------
# 1. Numerical Evaluation
# --------------------------------------------------

print("=" * 60)
print("GAUSS HYPERGEOMETRIC FUNCTION")
print("=" * 60)

a = 1.0
b = 2.0
c = 3.0

z_values = [0.0, 0.1, 0.2, 0.3, 0.5, 0.8]

print("\nNumerical evaluation of 2F1(a,b;c;z)")
print(f"Parameters: a={a}, b={b}, c={c}\n")

for z in z_values:
    value = hyp2f1(a, b, c, z)
    print(f"z = {z:4.2f}    2F1 = {value:.12f}")
# --------------------------------------------------
# 2. Analytical Verification
# --------------------------------------------------

print("\n" + "=" * 60)
print("ANALYTICAL VERIFICATION")
print("=" * 60)

def analytical_formula(z):
    if z == 0:
        return 1.0
    return 2 * (-np.log(1 - z) - z) / z**2


print("\nComparison between SciPy and analytical formula\n")

print(f"{'z':>8} {'SciPy':>18} {'Analytical':>18} {'Abs. Error':>15}")

for z in z_values:
    scipy_value = hyp2f1(a, b, c, z)
    exact_value = analytical_formula(z)
    error = abs(scipy_value - exact_value)

    print(
        f"{z:8.2f} "
        f"{scipy_value:18.12f} "
        f"{exact_value:18.12f} "
        f"{error:15.3e}"
    )

# --------------------------------------------------
# 3. Series Representation and Convergence
# --------------------------------------------------

print("\n" + "=" * 60)
print("SERIES REPRESENTATION AND CONVERGENCE")
print("=" * 60)

def hypergeometric_partial_sum(z, N):
    """
    Compute the partial sum of:
    2F1(1,2;3;z) = sum [2*z^n/(n+2)]
    """
    return sum(2 * z**n / (n + 2) for n in range(N))


z_test = 0.8
exact_value = hyp2f1(1, 2, 3, z_test)

N_values = [5, 10, 20, 50, 100, 200]

print(f"\nFixed parameter: z = {z_test}")
print(f"Reference value: {exact_value:.12f}\n")

print(f"{'N':>8} {'Partial Sum':>20} {'Absolute Error':>20}")

for N in N_values:
    approximation = hypergeometric_partial_sum(z_test, N)
    error = abs(approximation - exact_value)

    print(f"{N:8d} {approximation:20.12f} {error:20.3e}")
# --------------------------------------------------
# 4. Graphical Visualization
# --------------------------------------------------

print("\n" + "=" * 60)
print("GRAPHICAL VISUALIZATION")
print("=" * 60)

z_plot = np.linspace(0, 0.95, 300)

parameter_sets = [
    (1, 2, 3),
    (1, 1, 2),
    (2, 2, 4),
    (1, 3, 4)
]

plt.figure(figsize=(10, 6))

for a_plot, b_plot, c_plot in parameter_sets:
    values = hyp2f1(a_plot, b_plot, c_plot, z_plot)

    plt.plot(
        z_plot,
        values,
        label=f"a={a_plot}, b={b_plot}, c={c_plot}"
    )

plt.xlabel("z")
plt.ylabel(r"$\,_2F_1(a,b;c;z)$")
plt.title("Gauss Hypergeometric Function")
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig(
    "special_functions/hypergeometric_function/images/hypergeometric_curves.png",
    dpi=300
)

plt.show()

print("Graph saved successfully.")

# --------------------------------------------------
# 5. Gauss Hypergeometric Differential Equation
# --------------------------------------------------

print("\n" + "=" * 60)
print("HYPERGEOMETRIC DIFFERENTIAL EQUATION")
print("=" * 60)


def hypergeometric_derivative(a, b, c, z):
    return (a * b / c) * hyp2f1(a + 1, b + 1, c + 1, z)


def hypergeometric_second_derivative(a, b, c, z):
    coefficient = (a * (a + 1) * b * (b + 1)) / (c * (c + 1))
    return coefficient * hyp2f1(a + 2, b + 2, c + 2, z)


print("\nVerification of the differential equation")
print("z(1-z)y'' + [c-(a+b+1)z]y' - ab y = 0\n")

print(f"{'z':>8} {'Residual':>20}")

for z in [0.1, 0.2, 0.3, 0.5, 0.7, 0.9]:

    y = hyp2f1(a, b, c, z)
    y_prime = hypergeometric_derivative(a, b, c, z)
    y_double_prime = hypergeometric_second_derivative(a, b, c, z)

    residual = (
        z * (1 - z) * y_double_prime
        + (c - (a + b + 1) * z) * y_prime
        - a * b * y
    )

    print(f"{z:8.2f} {residual:20.3e}")
# --------------------------------------------------
# 6. Euler's Transformation Formula
# --------------------------------------------------

print("\n" + "=" * 60)
print("EULER'S TRANSFORMATION FORMULA")
print("=" * 60)

a_e, b_e, c_e = 1.5, 0.8, 2.5

print("\nVerification of Euler's transformation\n")
print(f"Parameters: a={a_e}, b={b_e}, c={c_e}")
print(f"{'z':>8} {'Left Side':>18} {'Right Side':>18} {'Abs. Error':>15}")

for z in [0.1, 0.2, 0.3, 0.5, 0.7]:

    left_side = hyp2f1(a_e, b_e, c_e, z)

    right_side = (
        (1 - z)**(c_e - a_e - b_e)
        * hyp2f1(c_e - a_e, c_e - b_e, c_e, z)
    )

    error = abs(left_side - right_side)

    print(
        f"{z:8.2f} "
        f"{left_side:18.12f} "
        f"{right_side:18.12f} "
        f"{error:15.3e}"
    )

# --------------------------------------------------
# 7. Euler's Integral Representation
# --------------------------------------------------

print("\n" + "=" * 60)
print("EULER'S INTEGRAL REPRESENTATION")
print("=" * 60)

a_i = 1.5
b_i = 1.2
c_i = 3.0

print(f"\nParameters: a={a_i}, b={b_i}, c={c_i}")
print(f"{'z':>8} {'SciPy':>18} {'Quadrature':>18} {'Abs. Error':>15}")

for z in [0.1, 0.2, 0.3, 0.5, 0.7]:

    def integrand(t):
        return (
            t**(b_i - 1)
            * (1 - t)**(c_i - b_i - 1)
            * (1 - z*t)**(-a_i)
        )

    beta_factor = math.gamma(c_i) / (
        math.gamma(b_i) * math.gamma(c_i - b_i)
    )

    integral_value, integration_error = quad(integrand, 0, 1)

    quadrature_value = beta_factor * integral_value
    scipy_value = hyp2f1(a_i, b_i, c_i, z)

    absolute_error = abs(scipy_value - quadrature_value)

    print(
        f"{z:8.2f} "
        f"{scipy_value:18.12f} "
        f"{quadrature_value:18.12f} "
        f"{absolute_error:15.3e}"
    )

# --------------------------------------------------
# 8. Special Values and Gauss's Summation Theorem
# --------------------------------------------------

print("\n" + "=" * 60)
print("SPECIAL VALUES AND GAUSS'S SUMMATION THEOREM")
print("=" * 60)

# Special value at z = 0
print("\n1. Special value at z = 0")

for parameters in [(1, 2, 3), (1.5, 0.8, 2.5), (2, 1, 4)]:
    a_s, b_s, c_s = parameters
    value = hyp2f1(a_s, b_s, c_s, 0)

    print(
        f"a={a_s}, b={b_s}, c={c_s}: "
        f"2F1(0) = {value:.12f}"
    )

# Gauss summation theorem at z = 1
print("\n2. Gauss summation theorem at z = 1")

parameter_sets = [
    (1.0, 1.0, 3.0),
    (0.5, 1.0, 2.5),
    (1.2, 0.8, 3.0),
    (1.5, 0.5, 3.0)
]

print(f"{'Parameters (a,b,c)':>25} {'SciPy':>18} "
      f"{'Gamma formula':>18} {'Abs. Error':>15}")

for a_s, b_s, c_s in parameter_sets:

    scipy_value = hyp2f1(a_s, b_s, c_s, 1.0)

    gamma_value = (
        math.gamma(c_s) * math.gamma(c_s - a_s - b_s)
        / (
            math.gamma(c_s - a_s)
            * math.gamma(c_s - b_s)
        )
    )

    error = abs(scipy_value - gamma_value)

    print(
        f"({a_s:.1f},{b_s:.1f},{c_s:.1f})"
        f"{scipy_value:18.12f}"
        f"{gamma_value:18.12f}"
        f"{error:15.3e}"
    )

# --------------------------------------------------
# 9. Accuracy and Convergence Near z = 1
# --------------------------------------------------

print("\n" + "=" * 60)
print("ACCURACY AND CONVERGENCE NEAR z = 1")
print("=" * 60)

z_values_near_one = [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]
N = 5000

print(f"\nNumber of terms: N = {N}\n")
print(f"{'z':>8} {'SciPy':>18} {'Partial Sum':>18} {'Abs. Error':>15}")

for z in z_values_near_one:

    scipy_value = hyp2f1(1, 2, 3, z)

    partial_sum = hypergeometric_partial_sum(z, N)

    error = abs(scipy_value - partial_sum)

    print(
        f"{z:8.2f} "
        f"{scipy_value:18.12f} "
        f"{partial_sum:18.12f} "
        f"{error:15.3e}"
    )

# --------------------------------------------------
# 10. Graphical Error Analysis
# --------------------------------------------------

print("\n" + "=" * 60)
print("GRAPHICAL ERROR ANALYSIS")
print("=" * 60)

z_error = np.linspace(0.1, 0.99, 150)

errors = []

for z in z_error:
    scipy_value = hyp2f1(1, 2, 3, z)
    partial_sum = hypergeometric_partial_sum(z, 100)

    errors.append(abs(scipy_value - partial_sum))

plt.figure(figsize=(10, 6))

plt.semilogy(z_error, errors, linewidth=2)

plt.xlabel("z")
plt.ylabel("Absolute Error")
plt.title("Convergence Error of the Hypergeometric Series")
plt.grid(True, which="both")
plt.tight_layout()

plt.savefig(
    "special_functions/hypergeometric_function/images/hypergeometric_convergence.png",
    dpi=300
)

plt.show()

print("Convergence error graph saved successfully.")
