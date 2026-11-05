
import numpy as np
import matplotlib.pyplot as plt

# Define the domain
x = np.linspace(-3, 3, 200)

# Define the exponential function
y = np.exp(x)

# Plot the function
plt.figure(figsize=(9, 5))
plt.plot(x, y, label=r'$f(x)=e^x$')

plt.title("Graph of the Exponential Function")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.legend()

plt.savefig("exponential_function.png")
plt.show()
