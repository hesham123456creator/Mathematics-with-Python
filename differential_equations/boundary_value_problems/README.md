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

## Numerical Methods and Tools

The projects use the following methods and libraries:

- Finite Difference Method.
- Tridiagonal linear system solvers.
- Numerical boundary value solvers.
- Nonlinear differential equation solvers.
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
- Scientific visualization using Matplotlib.

## Author

**Hesham Moustafa**

Assistant Professor of Mathematics

Fahd bin Sultan University  
Tabuk, Saudi Arabia

GitHub: https://github.com/hesham123456creator

Research interests: Numerical Analysis, Differential Equations, Computational Mathematics, and Scientific Computing.
