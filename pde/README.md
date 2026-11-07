# Partial Differential Equations Using Python

This directory contains Python implementations of numerical methods for solving partial differential equations (PDEs), with an emphasis on mathematical modeling, numerical stability, error analysis, and convergence.

## 1. Heat Diffusion Equation

The project considers the one-dimensional heat equation:

\[
\frac{\partial u}{\partial t}
=
\alpha \frac{\partial^2 u}{\partial x^2},
\qquad 0 < x < L,\quad t > 0.
\]

where:
- \(u(x,t)\) represents the temperature distribution.
- \(\alpha\) is the thermal diffusivity.
- \(x\) is the spatial variable.
- \(t\) is time.

### Initial and Boundary Conditions

The numerical example uses:

\[
u(x,0)=\sin(\pi x),
\]

with homogeneous Dirichlet boundary conditions:

\[
u(0,t)=u(1,t)=0.
\]

For \(L=1\), the exact analytical solution is:

\[
u(x,t)=e^{-\alpha\pi^2t}\sin(\pi x).
\]

## 2. Numerical Method

The heat equation is solved using the Forward-Time Central-Space (FTCS) finite difference method.

The discretization formula is:

\[
u_i^{n+1}
=
u_i^n+r(u_{i+1}^n-2u_i^n+u_{i-1}^n),
\]

where:

\[
r=\frac{\alpha\Delta t}{(\Delta x)^2}.
\]

For this explicit scheme, the stability condition is:

\[
0<r\leq\frac{1}{2}.
\]

## 3. Python Program

| File | Description |
|---|---|
| `heat_diffusion_project.py` | Numerical solution of the one-dimensional heat equation using FTCS, with analytical error verification and convergence analysis. |

## 4. Numerical Analysis Features

The implementation includes:

- Explicit finite difference discretization.
- Stability condition verification.
- Vectorized numerical computation using NumPy.
- Comparison with the exact analytical solution.
- Maximum, mean, and RMS error calculations.
- Temperature evolution visualization.
- Spatial grid refinement and convergence study.
- Observed convergence order calculation.
- High-resolution plot generation.

## 5. Required Libraries

- NumPy
- Matplotlib

Install the required libraries using:

```bash
pip install numpy matplotlib
```

Run the program using:

```bash
python heat_diffusion_project.py
```

## 6. Expected Output

The program generates three figures:

1. `heat_diffusion_evolution.png` – comparison of numerical and exact temperature distributions at selected times.
2. `heat_diffusion_error.png` – absolute numerical error at the final time.
3. `heat_diffusion_convergence.png` – log-log convergence plot showing how numerical error changes with grid spacing.

## 7. Learning Outcomes

This project demonstrates practical skills in:

- Numerical solutions of partial differential equations.
- Finite difference methods.
- Stability analysis of explicit numerical schemes.
- Error estimation and convergence analysis.
- Scientific programming and visualization using Python.

## Author

**Hesham Moustafa**  
Assistant Professor of Mathematics  
Fahd bin Sultan University, Saudi Arabia

**Research Interests:** Numerical Analysis, Differential Equations, Computational Mathematics, and Special Functions.
