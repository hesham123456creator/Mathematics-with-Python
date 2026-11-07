# Numerical Methods for Differential Equations Using Python

## Overview

This folder contains Python implementations of numerical methods for solving ordinary differential equations (ODEs). The projects demonstrate classical numerical algorithms, adaptive integration techniques, accuracy analysis, convergence studies, and graphical visualization.

The examples are designed to connect mathematical theory with practical scientific computing using Python.

## Projects

### 1. Euler's Method and Runge–Kutta Methods

| File | Description |
|---|---|
| `euler_vs_rk4.py` | Compares Euler's method and classical RK4 for a first-order ODE. |
| `euler_vs_rk4_comparison.py` | Compares Euler and RK4 using a problem with a known analytical solution, including error analysis. |

### 2. Classical Runge–Kutta Methods

| File | Description |
|---|---|
| `first_order_ode_rk4.py` | Solves a first-order ODE using classical fourth-order Runge–Kutta. |
| `runge_kutta_rk4.py` | Implements RK4 with accuracy and convergence analysis. |
| `second_order_ode_rk4.py` | Solves a homogeneous second-order ODE by converting it into a first-order system. |
| `nonhomogeneous_ode_rk4.py` | Solves a nonhomogeneous second-order ODE using RK4. |
| `systems_of_odes_rk4.py` | Implements RK4 for a system of first-order differential equations. |

### 3. Adaptive Runge–Kutta Methods

| File | Description |
|---|---|
| `runge_kutta_rk45.py` | Solves a first-order ODE using SciPy's adaptive RK45 method. |
| `runge_kutta_rk45_example2.py` | Applies adaptive RK45 to an ODE with a known analytical solution and studies the numerical error. |
| `systems_of_odes.py` | Solves a system of ODEs using SciPy's adaptive RK45 solver. |

## Numerical Methods and Techniques

The projects cover:

- Euler's numerical integration method.
- Classical fourth-order Runge–Kutta (RK4).
- Adaptive Runge–Kutta method (RK45).
- Conversion of higher-order ODEs into first-order systems.
- Numerical solution of systems of differential equations.
- Absolute error and maximum error analysis.
- Step-size selection and adaptive integration.
- Fourth-order convergence studies.
- Phase-plane visualization.
- Energy conservation analysis.

## Python Libraries

The following libraries are used:

- **NumPy:** Numerical arrays and mathematical computations.
- **Matplotlib:** Graphical visualization and plotting.
- **SciPy:** Adaptive numerical integration using `solve_ivp`.

## Installation

Install the required libraries:

```bash
pip install numpy scipy matplotlib
```

Alternatively, install all repository dependencies:

```bash
pip install -r requirements.txt
```

## Running the Programs

Navigate to this folder and execute a Python file:

```bash
python runge_kutta_rk4.py
```

For example, to run the adaptive RK45 solver:

```bash
python runge_kutta_rk45.py
```

Each program prints numerical results and error information and generates relevant plots when executed.

## Learning Outcomes

Through these projects, users can:

1. Understand the implementation of numerical integration methods.
2. Compare the accuracy of Euler's method and Runge–Kutta methods.
3. Solve first-order and higher-order ordinary differential equations numerically.
4. Apply numerical methods to systems of differential equations.
5. Investigate numerical stability and convergence behavior.
6. Evaluate numerical solutions against exact analytical solutions.
7. Visualize numerical solutions, errors, and phase-plane trajectories.
8. Use scientific Python libraries to solve mathematical problems.

## Author

**Hesham Moustafa**

Assistant Professor of Mathematics  
Fahd bin Sultan University, Tabuk, Saudi Arabia

Research interests include numerical analysis, differential equations, computational mathematics, and special functions.

GitHub: [hesham123456creator](https://github.com/hesham123456creator)

Repository: [Mathematics-with-Python](https://github.com/hesham123456creator/Mathematics-with-Python)
