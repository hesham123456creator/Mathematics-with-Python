# Gamma Function: Mathematical Properties and Computation

## 1. Overview

This project explores the Gamma function using Python. It combines mathematical properties, numerical computation, verification of important identities, and graphical visualization.

The Gamma function is a fundamental special function that extends the factorial function to non-integer arguments.

## 2. Mathematical Definition

For \(x>0\), the Gamma function is defined by

\[
\Gamma(x)=\int_0^\infty t^{x-1}e^{-t}\,dt.
\]

## 3. Fundamental Properties

### 3.1 Recurrence Relation

\[
\Gamma(x+1)=x\Gamma(x).
\]

### 3.2 Factorial Relation

For positive integers \(n\),

\[
\Gamma(n)=(n-1)!.
\]

### 3.3 Special Value

\[
\Gamma\left(\frac12\right)=\sqrt{\pi}.
\]

### 3.4 Minimum on the Positive Real Axis

The Gamma function has a unique global minimum for positive real arguments, approximately at

\[
x_{\min}\approx 1.461632.
\]

The corresponding minimum value is

\[
\Gamma(x_{\min})\approx 0.885603.
\]

## 4. Python Libraries

The implementation uses:

- NumPy
- SciPy
- SymPy
- Matplotlib
- Math (Python standard library)

## 5. Numerical Experiments

The Python program performs:

1. Numerical evaluation of the Gamma function.
2. Verification of the recurrence relation.
3. Verification of the factorial relation.
4. Verification of the special value at \(x=1/2\).
5. Graphical visualization of the Gamma function.
6. Numerical determination of its minimum using SciPy optimization.

## 6. Visualizations

### 6.1 Gamma Function

![Gamma Function](gamma_function_plot.png)

### 6.2 Gamma Function and Its Numerical Minimum

![Gamma Function Minimum](gamma_function_minimum.png)

## 7. Project Structure

```text
gamma_function/
├── gamma_function.py
├── gamma_function_plot.png
├── gamma_function_minimum.png
└── README.md

## 8. How to Run

From the root directory of the repository:

```bash
python special_functions/gamma_function/gamma_function.py
```

## 9. Learning Outcomes

This project demonstrates practical skills in:

- Special functions and their mathematical properties.
- Numerical evaluation using SciPy.
- Verification of mathematical identities.
- Numerical optimization.
- Scientific computing with Python.
- Mathematical visualization using Matplotlib.

## Author

**Hesham Moustafa**

Assistant Professor of Mathematics
