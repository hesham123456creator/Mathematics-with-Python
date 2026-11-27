import mpmath as mp
import matplotlib.pyplot as plt
import os

# --------------------------------------------------
# Rigorous bounds for the Gauss-Psi remainder
#
# R_N(z) = sum_{n=N+1}^\infty a_n z^n
# a_n = 2 H_n / (n+2)
#
# L_p(N,z) < R_N(z) < U_p(N,z)
# --------------------------------------------------

mp.mp.dps = 100

os.makedirs("images", exist_ok=True)


def coefficient(n):
    """Compute a_n = 2 H_n / (n+2)."""
    return 2 * mp.harmonic(n) / (n + 2)


def exact_function(z):
    """Closed-form representation of D(z)."""
    return (2 / z**2) * (
        -z
        + mp.log(1 - z)**2 / 2
        - (1 - z) * mp.log(1 - z)
    )


def partial_sum(N, z):
    """Compute the first N terms of D(z)."""
    return mp.fsum(
        coefficient(n) * z**n
        for n in range(1, N + 1)
    )


def actual_remainder(N, z):
    """High-precision evaluation of R_N(z)."""
    return exact_function(z) - partial_sum(N, z)


def lower_bound(N, z, p):
    """Lower bound L_p(N,z)."""
    return mp.fsum(
        coefficient(N + 1 + k) * z**(N + 1 + k)
        for k in range(p)
    )


def upper_bound(N, z, p):
    """Upper bound U_p(N,z)."""
    L = lower_bound(N, z, p)

    return (
        L
        + coefficient(N + p + 1)
        * z**(N + p + 1) / (1 - z)
    )


# --------------------------------------------------
# Numerical verification
# --------------------------------------------------

print("\nGAUSS-PSI REMAINDER BOUNDS")
print("=" * 75)

N_values = [10, 50, 100]
z_values = [0.3, 0.5, 0.7, 0.9]
p_values = [1, 2, 5, 10, 20]

for z in z_values:
    print(f"\n z = {z}")

    for N in N_values:
        R = actual_remainder(N, mp.mpf(str(z)))

        print(f"\n N = {N}")
        print(f"Actual remainder = {mp.nstr(R, 12)}")

        for p in p_values:
            L = lower_bound(N, mp.mpf(str(z)), p)
            U = upper_bound(N, mp.mpf(str(z)), p)

            verified = L < R < U

            width = U - L

            print(
                f"p={p:2d} | "
                f"L={mp.nstr(L, 8):>14} | "
                f"U={mp.nstr(U, 8):>14} | "
                f"Width={mp.nstr(width, 5):>12} | "
                f"Verified={verified}"
            )


# --------------------------------------------------
# Plot: Enclosure width versus p
# --------------------------------------------------

N = 50
z = mp.mpf("0.7")

p_plot = list(range(1, 31))
widths = [
    float(upper_bound(N, z, p) - lower_bound(N, z, p))
    for p in p_plot
]

plt.figure(figsize=(9, 6))

plt.semilogy(p_plot, widths, marker="o")

plt.xlabel("Number of retained terms p")
plt.ylabel("Width of the enclosure")
plt.title("Convergence of Gauss-Psi Remainder Bounds")

plt.grid(True, which="both", linestyle="--", alpha=0.6)
plt.tight_layout()

plt.savefig(
    "images/gauss_psi_remainder_bounds.png",
    dpi=300
)

plt.show()

print("\nPlot saved successfully!")
