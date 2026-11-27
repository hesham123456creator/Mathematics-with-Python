"""
Parameter Study of the Lerch Transcendent
=========================================
Investigates the influence of z, s, and a.
"""

import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp

mp.mp.dps = 30


def calculate_lerch(z, s, a):
    return float(mp.lerchphi(z, s, a))


# --------------------------------------------------
# Experiment 1: Vary z
# --------------------------------------------------

z_values = np.linspace(-0.9, 0.9, 100)
s_fixed = 2
a_fixed = 1

z_results = [
    calculate_lerch(z, s_fixed, a_fixed)
    for z in z_values
]


# --------------------------------------------------
# Experiment 2: Vary s
# --------------------------------------------------

s_values = np.linspace(0.5, 5, 100)
z_fixed = 0.5
a_fixed = 1

s_results = [
    calculate_lerch(z_fixed, s, a_fixed)
    for s in s_values
]


# --------------------------------------------------
# Experiment 3: Vary a
# --------------------------------------------------

a_values = np.linspace(0.2, 5, 100)
z_fixed = 0.5
s_fixed = 2

a_results = [
    calculate_lerch(z_fixed, s_fixed, a)
    for a in a_values
]


# --------------------------------------------------
# Plot 1: Dependence on z
# --------------------------------------------------

plt.figure(figsize=(9, 6))
plt.plot(z_values, z_results)
plt.xlabel("z")
plt.ylabel("Lerch Phi(z, 2, 1)")
plt.title("Lerch Transcendent: Variation with z")
plt.grid(True, alpha=0.4)
plt.tight_layout()

plt.savefig(
    "special_functions/lerch_function/images/"
    "lerch_variation_z.png",
    dpi=300
)


# --------------------------------------------------
# Plot 2: Dependence on s
# --------------------------------------------------

plt.figure(figsize=(9, 6))
plt.plot(s_values, s_results)
plt.xlabel("s")
plt.ylabel("Lerch Phi(0.5, s, 1)")
plt.title("Lerch Transcendent: Variation with s")
plt.grid(True, alpha=0.4)
plt.tight_layout()

plt.savefig(
    "special_functions/lerch_function/images/"
    "lerch_variation_s.png",
    dpi=300
)


# --------------------------------------------------
# Plot 3: Dependence on a
# --------------------------------------------------

plt.figure(figsize=(9, 6))
plt.plot(a_values, a_results)
plt.xlabel("a")
plt.ylabel("Lerch Phi(0.5, 2, a)")
plt.title("Lerch Transcendent: Variation with a")
plt.grid(True, alpha=0.4)
plt.tight_layout()

plt.savefig(
    "special_functions/lerch_function/images/"
    "lerch_variation_a.png",
    dpi=300
)

plt.show()

print("Parameter study completed successfully.")
print("Three graphical results have been generated.")
