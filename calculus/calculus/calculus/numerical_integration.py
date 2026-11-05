
import numpy as np
from scipy.integrate import quad

# Define the function
def f(x):
    return np.exp(-x**2)

# Calculate the definite integral
result, error = quad(f, 0, 1)

# Display results
print("Numerical integral =", result)
print("Estimated error =", error)
