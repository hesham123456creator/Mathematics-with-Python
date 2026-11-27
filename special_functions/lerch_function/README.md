# Lerch Transcendent

## 1. Introduction

The Lerch transcendent is an important special function that generalizes several classical mathematical functions.

It is defined by

\[
\Phi(z,s,a)=
\sum_{n=0}^{\infty}\frac{z^n}{(n+a)^s}.
\]

For the numerical experiments in this project, we consider real parameters satisfying

\[
|z|<1,\qquad s>0,\qquad a>0.
\]

## 2. Mathematical Background

The Lerch transcendent is closely related to several important special functions, including:

- Hurwitz zeta function.
- Polylogarithm function.
- Riemann zeta function.

Its integral representation is

\[
\Phi(z,s,a)=
\frac{1}{\Gamma(s)}
\int_0^\infty
\frac{t^{s-1}e^{-at}}
{1-ze^{-t}}\,dt.
\]

## 3. Numerical Implementations

Two independent numerical approaches are implemented.

### 3.1 Series Representation

File: `lerch_function.py`

The infinite series is approximated using a finite number of terms:

\[
\Phi(z,s,a)\approx
\sum_{n=0}^{N-1}
\frac{z^n}{(n+a)^s}.
\]

### 3.2 Integral Representation

File: `integral_representation.py`

The integral is evaluated numerically using SciPy's adaptive quadrature algorithm.

## 4. Numerical Verification

File: `numerical_verification.py`

The numerical implementation is compared with the independent reference function `mpmath.lerchphi()`.

For the test parameters

\[
z=0.5,\quad s=2,\quad a=1,
\]

the reference value is

\[
\Phi(0.5,2,1)
\approx 1.164481052930025.
\]

Both numerical implementations agree with the reference to double-precision accuracy.

## 5. Convergence Analysis

The convergence behavior is investigated using different truncation levels.

The experiments distinguish between:

- Truncation error.
- Floating-point error.

The results demonstrate rapid convergence of the series and the limitations of double-precision arithmetic.

Figure:

`images/lerch_convergence.png`

## 6. Parameter Study

File: `parameter_study.py`

Three numerical experiments investigate the influence of the parameters.

### Variation with z

`images/lerch_variation_z.png`

### Variation with s

`images/lerch_variation_s.png`

### Variation with a

`images/lerch_variation_a.png`

## 7. Mathematical Properties

File: `property_verification.py`

### 7.1 Monotonicity with respect to z

For real parameters satisfying

\[
0<z<1,\qquad s>0,\qquad a>0,
\]

the Lerch transcendent is strictly increasing with respect to z:

\[
\frac{\partial\Phi(z,s,a)}{\partial z}
=
\sum_{n=1}^{\infty}
\frac{n z^{n-1}}{(n+a)^s}>0.
\]

### 7.2 Monotonicity with respect to s

For

\[
0<z<1,\qquad a\geq1,\qquad s>0,
\]

the Lerch transcendent is strictly decreasing with respect to s:

\[
\frac{\partial\Phi(z,s,a)}{\partial s}
=
-\sum_{n=0}^{\infty}
\frac{z^n\ln(n+a)}{(n+a)^s}<0.
\]

### 7.3 Derivative with respect to a

An important analytical identity is

\[
\boxed{
\frac{\partial\Phi(z,s,a)}{\partial a}
=
-s\Phi(z,s+1,a)
}
\]

This identity is verified numerically using high-precision arithmetic and finite differences.

## 8. Project Structure

```text
lerch_function/
│
├── README.md
├── lerch_function.py
├── numerical_verification.py
├── integral_representation.py
├── parameter_study.py
├── property_verification.py
│
└── images/
    ├── lerch_convergence.png
    ├── lerch_variation_z.png
    ├── lerch_variation_s.png
    └── lerch_variation_a.png


```

## Technologies

The project uses the following Python libraries:

- **NumPy:** Numerical computations and array operations.
- **SciPy:** Numerical evaluation and scientific computing.
- **Matplotlib:** Visualization and graphical representation.
- **mpmath:** High-precision evaluation of special functions.

## Objectives

The main objectives of this project are:

1. Implement numerical evaluation of the Lerch transcendent.
2. Study its convergence properties.
3. Investigate the effects of the parameters \(z\), \(s\), and \(a\).
4. Compare numerical approximations with high-precision reference values.
5. Study the behavior of the function through numerical experiments.
6. Develop graphical representations of the Lerch transcendent.
7. Build a collection of Python tools for numerical investigations of special functions.


