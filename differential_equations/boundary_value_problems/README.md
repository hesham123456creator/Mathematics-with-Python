# Boundary Value Problems in Python

## Overview

This project explores numerical methods for solving linear and nonlinear boundary value problems (BVPs) using Python.

The main objectives are to implement finite difference schemes, use SciPy's numerical boundary value solver, compare numerical solutions with exact analytical solutions, and investigate numerical accuracy and convergence.

## Mathematical Background

A general second-order boundary value problem can be written as

\[
y''(x)=f(x,y,y'), \qquad a\leq x\leq b,
\]

subject to boundary conditions such as

\[
y(a)=\alpha,\qquad y(b)=\beta.
\]

This repository investigates several examples of such problems.

## Python Projects

### 1. `finite_difference_method.py`

Implements the Finite Difference Method for the boundary value problem

\[
y''(x)=-\pi^2\sin(\pi x),
\]

with boundary conditions

\[
y(0)=0,\qquad y(1)=0.
\]

The numerical solution is compared with the exact solution \(y(x)=\sin(\pi x)\).

### 2. `finite_difference_accuracy.py`

Investigates the accuracy of the Finite Difference Method using different grid sizes.

Features:
- Maximum absolute error calculations.
- Observed convergence order.
- Log-log convergence plots.
- Comparison with the theoretical second-order convergence rate.

### 3. `compare_bvp_methods.py`

Compares two numerical approaches for solving a linear boundary value problem:

- Finite Difference Method.
- SciPy `solve_bvp()`.

The numerical results are compared with the exact analytical solution.

### 4. `solve_bvp_example.py`

Demonstrates how to solve a second-order boundary value problem using SciPy's `solve_bvp()` function.

The project includes:
- Conversion into a first-order system.
- Boundary condition implementation.
- Numerical solution evaluation.
- Error analysis and boundary residual verification.

### 5. `nonlinear_bvp.py`

Introduces nonlinear boundary value problems through

\[
y''+y^2=0,
\]

with homogeneous boundary conditions.

The example demonstrates the trivial solution and numerical verification using SciPy.

### 6. `nonlinear_bvp_nontrivial.py`

Studies the nonlinear Bratu boundary value problem:

\[
y''+\lambda e^y=0,
\]

subject to

\[
y(0)=0,\qquad y(1)=0.
\]

The numerical solution is compared with an exact analytical solution.

This example includes numerical error analysis, boundary condition verification, and differential equation residual checking.

### 7. `bratu_parameter_study.py`

Investigates both analytical solution branches of the nonlinear Bratu boundary value problem:

\[
y''+\lambda e^y=0,
\]

subject to

\[
y(0)=0,\qquad y(1)=0.
\]

The study considers four parameter values:

\[
\lambda=1.0,\quad 2.0,\quad 3.0,\quad 3.5.
\]

For each parameter value below the critical value, two distinct positive solutions exist: a lower branch and an upper branch.

The critical parameter is approximately

\[
\lambda_c=3.5138307191.
\]

The numerical solutions obtained using SciPy's `solve_bvp()` are compared with the exact analytical solutions.

Features:
- Computation of both analytical solution branches.
- Numerical branch tracking using suitable initial guesses.
- Critical parameter estimation using `brentq()`.
- Maximum absolute error analysis.
- Differential equation residual evaluation.
- Collocation RMS residual analysis.
- Comparison of both branches through scientific visualization.

The resulting plot is saved as:

`images/bratu_parameter_study.png`

The numerical results demonstrate excellent accuracy, with maximum absolute errors below \(1.1\times10^{-10}\).
Investigates both analytical solution branches of the nonlinear Bratu boundary value problem:

\[
y''+\lambda e^y=0,
\]

subject to

\[
y(0)=0,\qquad y(1)=0.
\]

The study considers four parameter values:

\[
\lambda=1.0,\quad 2.0,\quad 3.0,\quad 3.5.
\]

For each parameter value below the critical value, two distinct positive solutions exist: a lower branch and an upper branch.

The critical parameter is approximately

\[
\lambda_c=3.5138307191.
\]

The numerical solutions obtained using SciPy's `solve_bvp()` are compared with the exact analytical solutions.

Features:
- Computation of both analytical solution branches.
- Numerical branch tracking using suitable initial guesses.
- Critical parameter estimation using `brentq()`.
- Maximum absolute error analysis.
- Differential equation residual evaluation.
- Collocation RMS residual analysis.
- Comparison of both branches through scientific visualization.

The resulting plot is saved as:

`images/bratu_parameter_study.png`

The numerical results demonstrate excellent accuracy, with maximum absolute errors below \(1.1\times10^{-10}\).

### 8. `nonlinear_bvp_parameter.py`

Demonstrates the numerical solution of a nonlinear boundary value problem involving an unknown parameter:

\[
y''+\lambda y^2=0,
\]

subject to

\[
y(0)=0,\qquad y'(0)=1,\qquad y(1)=0.5.
\]

The unknown parameter is estimated using SciPy's `solve_bvp()` function.

The computed parameter is approximately

\[
\lambda=8.551517046259.
\]

Features:
- Simultaneous estimation of an unknown parameter and numerical solution.
- Adaptive mesh refinement.
- Collocation residual analysis.
- Differential equation residual evaluation.
- Boundary condition verification.
- Numerical solution visualization.

## Numerical Methods and Tools

The projects use the following methods and libraries:

- Finite Difference Method.
- Tridiagonal linear system solvers.
- Numerical boundary value solvers.
- Nonlinear differential equation solvers.
- Parameter estimation in nonlinear BVPs.
- Convergence and error analysis.
- NumPy.
- SciPy.
- Matplotlib.

## Installation

Install the required Python libraries:

```bash
pip install numpy scipy matplotlib
```

## Running the Programs

Run any Python program from this directory using:

```bash
python finite_difference_method.py
```

Replace the filename with the program you want to execute.

## Learning Outcomes

Through these projects, the repository demonstrates:

- Numerical solutions of linear and nonlinear BVPs.
- Implementation of finite difference discretization.
- Use of SciPy for boundary value problems.
- Comparison between numerical and analytical solutions.
- Numerical accuracy and convergence analysis.
- Estimation of unknown parameters in nonlinear differential equations.
- Scientific visualization using Matplotlib.

## Author

**Hesham Moustafa**

Assistant Professor of Mathematics

Fahd bin Sultan University  
Tabuk, Saudi Arabia

GitHub: https://github.com/hesham123456creator

Research interests: Numerical Analysis, Differential Equations, Computational Mathematics, and Scientific Computing.
