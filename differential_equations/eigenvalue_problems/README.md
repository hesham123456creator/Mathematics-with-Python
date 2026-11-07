# Numerical Eigenvalue Problems

## Overview

This project investigates eigenvalue problems associated with second-order differential equations using numerical methods in Python.

The main objective is to compare numerical approximations with exact analytical solutions and study the accuracy of finite difference methods.

## Mathematical Problem

Consider the boundary value problem:

\[
-y''(x)=\lambda y(x), \quad 0<x<1
\]

subject to:

\[
y(0)=0,\qquad y(1)=0.
\]

The exact eigenvalues are:

\[
\lambda_n=n^2\pi^2,\qquad n=1,2,3,\ldots
\]

The corresponding eigenfunctions are:

\[
y_n(x)=\sin(n\pi x).
\]

## Python Programs

### 1. eigenvalue_problem.py

Computes numerical approximations of eigenvalues using finite difference discretization.

### 2. eigenvalue_accuracy.py

Investigates the accuracy of numerical eigenvalues compared with exact analytical values.

### 3. compare_eigenfunctions.py

Compares numerical eigenfunctions with exact eigenfunctions and calculates maximum absolute errors.

## Numerical Methods and Libraries

- Finite Difference Method
- NumPy
- SciPy
- Matplotlib

## Expected Results

The numerical eigenvalues and eigenfunctions should converge toward their exact analytical counterparts as the mesh is refined.

The project demonstrates the relationship between differential equations, linear algebra, and numerical analysis.

## Author

Hesham Moustafa

Assistant Professor of Mathematics

GitHub: https://github.com/hesham123456creator
