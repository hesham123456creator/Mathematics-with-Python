# Gauss–Psi Connection

## 1. Introduction

This project studies the parameter derivative of the Gauss hypergeometric function and its connection with the digamma function and harmonic numbers.

We consider

\[
D(z)=\left.\frac{\partial}{\partial a}
{}_2F_1(a,2;3;z)\right|_{a=1}.
\]

## 2. Harmonic-Series Representation

Using the parameter derivative of the Pochhammer symbol, we obtain

\[
D(z)=\sum_{n=1}^{\infty}
\frac{2H_n}{n+2}z^n,
\]

where \(H_n\) denotes the nth harmonic number.

## 3. Closed-Form Representation

The function also admits the following logarithmic representation:

\[
D(z)=\frac{2}{z^2}
\left[
-z+\frac12\ln^2(1-z)
-(1-z)\ln(1-z)
\right].
\]

For \(0<z<1\).

## 4. Asymptotic Remainder Analysis

Define the truncation remainder by

\[
R_N(z)=\sum_{n=N+1}^{\infty}
\frac{2H_n}{n+2}z^n.
\]

The leading-order asymptotic approximation is

\[
R_N(z)\sim
\frac{2H_{N+1}}{N+3}
\frac{z^{N+1}}{1-z}.
\]

A second-order approximation is also investigated numerically.

## 5. Numerical Implementation

The Python implementation compares:

- Numerical differentiation of the hypergeometric function.
- The harmonic-number series.
- The analytical closed form.
- Leading-order and second-order remainder approximations.

## 6. Numerical Results

### Comparison of Three Representations

![Comparison](images/gauss_psi_comparison.png)

### Truncation Error

![Truncation Error](images/gauss_psi_truncation_error.png)

### Asymptotic Ratio

![Asymptotic Ratio](images/gauss_psi_asymptotic_ratio.png)

### Asymptotic Accuracy

![Asymptotic Accuracy](images/gauss_psi_asymptotic_accuracy.png)

## 7. Technologies

- Python
- NumPy
- SciPy
- Matplotlib
- mpmath

## 8. Research Direction

Future work will focus on rigorous asymptotic estimates, higher-order remainder approximations, and analytical inequalities.
