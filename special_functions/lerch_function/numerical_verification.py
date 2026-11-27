"""
Numerical Verification of the Lerch Transcendent
================================================
Comparison with mpmath and convergence analysis.
"""

import mpmath as mp
import matplotlib.pyplot as plt

from lerch_function import lerch_phi


mp.mp.dps = 100

z = 0.5
s = 2
a = 1

reference = mp.lerchphi(z, s, a)

N_values = [5, 10, 20, 30, 50, 75, 100, 150, 200]

truncation_errors = []
floating_point_errors = []

for N in N_values:

    # High-precision truncated series
    partial_sum = mp.fsum(
        mp.mpf(z)**n / (mp.mpf(n + a)**s)
        for n in range(N)
    )

    truncation_error = abs(reference - partial_sum)

    # Standard double-precision implementation
    numerical_value = lerch_phi(z, s, a, N)

    floating_error = abs(
        float(numerical_value) - float(reference)
    )

    truncation_errors.append(float(truncation_error))
    floating_point_errors.append(floating_error)

    print(
        f"N = {N:4d} | "
        f"Truncation error = {float(truncation_error):.6e} | "
        f"Floating-point error = {floating_error:.6e}"
    )


# Plot convergence
plt.figure(figsize=(9, 6))

plt.loglog(
    N_values,
    truncation_errors,
    'o-',
    label='Truncation error'
)

plt.loglog(
    N_values,
    floating_point_errors,
    's--',
    label='Floating-point error'
)

plt.xlabel("Number of terms (N)")
plt.ylabel("Absolute error")
plt.title("Convergence Analysis of the Lerch Transcendent")

plt.grid(True, which="both", linestyle="--", alpha=0.6)
plt.legend()
plt.tight_layout()

plt.savefig(
    "special_functions/lerch_function/images/"
    "lerch_convergence.png",
    dpi=300
)

plt.show()
