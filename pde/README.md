# Partial Differential Equations Using Python

This directory contains Python implementations of numerical methods for solving partial differential equations (PDEs), with an emphasis on mathematical modeling, numerical stability, iterative methods, error analysis, and convergence.

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

## 3. Additional PDE Examples

The directory also includes independent programs illustrating other PDE models and numerical approaches.

| File | Description |
|---|---|
| `heat_equation.py` | Numerical solution of the one-dimensional heat equation using an explicit finite difference scheme. |
| `laplace_equation.py` | Numerical solution of the two-dimensional Laplace equation using the Gauss-Seidel iterative method. |
| `wave_equation.py` | Numerical solution of the one-dimensional wave equation using an explicit finite difference method. |
| `compare_iterative_methods.py` | Comparison of Jacobi and Gauss-Seidel methods for solving a two-dimensional Laplace problem. |

### Jacobi versus Gauss-Seidel

The comparison program investigates two iterative methods for solving the discretized Laplace equation. It examines:

- Number of iterations required for convergence.
- Execution time.
- Difference between the computed numerical solutions.

The example illustrates that fewer iterations do not necessarily imply shorter execution time, since implementation and vectorization also affect computational performance.

## 4. Advanced Heat Diffusion Project

| File | Description |
|---|---|
| `heat_diffusion_project.py` | Numerical solution of the one-dimensional heat equation using FTCS, with analytical error verification and convergence analysis. |

The advanced implementation includes:

- Stability condition verification.
- Vectorized numerical computation using NumPy.
- Comparison with the exact analytical solution.
- Maximum, mean, and RMS error calculations.
- Temperature evolution visualization.
- Spatial grid refinement and convergence study.
- Observed convergence order calculation.
- High-resolution plot generation.

## 5. Wave Equation

The project also considers the one-dimensional wave equation:

\[
\frac{\partial^2 u}{\partial t^2}
=
c^2\frac{\partial^2 u}{\partial x^2},
\qquad 0 < x < L,\quad t > 0.
\]

where:
- \(u(x,t)\) represents the wave displacement.
- \(c\) is the wave propagation speed.
- \(x\) is the spatial variable.
- \(t\) is time.

### Numerical Method

The wave equation is solved using an explicit central finite difference scheme:

\[
u_i^{n+1}
=
2u_i^n-u_i^{n-1}
+r^2(u_{i+1}^n-2u_i^n+u_{i-1}^n),
\]

where:

\[
r=\frac{c\Delta t}{\Delta x}.
\]

The stability condition for this scheme is:

\[
0<r\leq 1.
\]

### Accuracy and Convergence Analysis

The implementation includes:
- Stability parameter verification.
- Numerical accuracy analysis.
- Maximum absolute error calculation.
- Root mean square error calculation.
- Grid refinement study.
- Observed convergence order calculation.
- Visualization of numerical and exact solutions.

The numerical experiment gives a stability parameter of approximately \(0.9091\).

The grid refinement study demonstrates second-order convergence, with an observed order approaching \(2\).

### Generated Figures

The wave equation program generates:

1. `wave_equation_comparison.png` – comparison of numerical and exact wave solutions.
2. `wave_equation_convergence.png` – convergence behavior under successive grid refinements.

## 6. Required Libraries

- NumPy
- SciPy
- Matplotlib

Install the required libraries using:

```bash
pip install numpy scipy matplotlib
```

## 7. Running the Programs

From the repository root, run the individual programs using:

```bash
python pde/heat_equation.py
python pde/laplace_equation.py
python pde/wave_equation.py
python pde/compare_iterative_methods.py
python pde/heat_diffusion_project.py
```

## 8. Expected Output

The programs generate numerical results and visualizations illustrating:

- Heat distribution over time.
- Numerical solutions of Laplace's equation.
- Wave propagation.
- Convergence behavior of iterative methods.
- Numerical errors and convergence rates for the heat diffusion problem.

The advanced heat diffusion project generates three figures:

1. `heat_diffusion_evolution.png` – comparison of numerical and exact temperature distributions at selected times.
2. `heat_diffusion_error.png` – absolute numerical error at the final time.
3. `heat_diffusion_convergence.png` – log-log convergence plot showing how numerical error changes with grid spacing.

The wave equation program generates two figures:

1. `wave_equation_comparison.png` – comparison of numerical and exact wave solutions.
2. `wave_equation_convergence.png` – convergence plot showing the numerical error as the spatial grid is refined.

## 9. Learning Outcomes

This project demonstrates practical skills in:

- Numerical solutions of partial differential equations.
- Finite difference methods.
- Stability analysis of explicit numerical schemes.
- Iterative methods for elliptic PDEs.
- Comparison of numerical algorithms.
- Error estimation and convergence analysis.
- Scientific programming and visualization using Python.

## Author

**Hesham Moustafa**  
Assistant Professor of Mathematics  
Fahd bin Sultan University, Saudi Arabia

**Research Interests:** Numerical Analysis, Differential Equations, Computational Mathematics, and Special Functions.
