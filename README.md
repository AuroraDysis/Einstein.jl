# Einstein

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AuroraDysis.github.io/Einstein.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AuroraDysis.github.io/Einstein.jl/dev/)
[![Build Status](https://github.com/AuroraDysis/Einstein.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AuroraDysis/Einstein.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AuroraDysis/Einstein.jl/graph/badge.svg?token=C99DVUUULL)](https://codecov.io/gh/AuroraDysis/Einstein.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

> [!WARNING]  
> I'm still actively developing it and migrating code from private repositories to this one, so it may currently lack some features.

## Motivation and Purpose

After developing hundreds of PDE solvers for various projects, I decided to create `Einstein.jl` because I couldn't find a library that was both highly convenient and effective for solving partial differential equations (PDEs) in the context of general relativity (GR).

When solving GR PDEs in spherical symmetry or even axisymmetry, the requirements are often far more demanding than in 3+1 cases. These include the need for very long-time simulations, conservation of quantities, highly accurate late-time behavior, and precise computation of quasinormal modes. While I have addressed these challenges individually in the past, I grew tired of copying and pasting solutions from old repositories. Instead, I decided to develop a comprehensive, feature-rich, and highly efficient library. This is my answer to those challenges.

## Description

**Einstein.jl** is a **high-performance** suite designed to compute **arbitrary-precision solutions** of PDEs in GR. It is planned to support solving the following built-in equations (some of which were developed for previous projects and will be ported to this library):

### Spatial Discretization

- [x] Chebyshev collocation at Chebyshev points of the first and second kinds (Most of algorithms are translated from Chebfun)
- [x] Finite difference method
- [x] Hermite finite difference method
- [x] Rectangular collocation method (Most of algorithms are translated from Chebfun)
- [x] Ultraspherical spectral method (Most of algorithms are translated from Chebfun)
  - For boundary value problems, I recommend using [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl).
- [ ] Ultraspherical rectangular collocation ([tomtrogdon/URCMethod.jl: A Jupyter notebook with an implementation of the URC method](https://github.com/tomtrogdon/URCMethod.jl))

### General Utilities

- [x] Correctly rounded floating-point dot/sum using [xsum](https://arxiv.org/abs/1505.05571) or [Kahan summation](https://en.wikipedia.org/wiki/Kahan_summation_algorithm)
- [x] Spin-weighted spheroidal harmonics using the Cook-Zalutskiy spectral method

### Quasinormal Modes

- [ ] Compute quasinormal modes for the Schwarzschild black hole using the Regge-Wheeler-Zerilli equation.
  - [ ] Continued fraction method to determine the eigenvalue
  - [x] Ultraspherical spectral method in hyperboloidal coordinates to determine the eigenfunction.
  - [x] Dolan and Ottewill Regge poles expansion method
  - [ ] WKB approximation
  - [ ] Direct integration method
- [x] Compute quasinormal modes for the Kerr black hole using the Teukolsky equation.
  - [x] Continued fraction method for the radial equation (translated from [duetosymmetry/qnm](https://github.com/duetosymmetry/qnm))
  - [x] Cook-Zalutskiy spectral method for the angular sector (translated from [duetosymmetry/qnm](https://github.com/duetosymmetry/qnm))
  - [x] Direct integration method for the radial equation
  - [x] Ultraspherical spectral method in hyperboloidal coordinates to determine the eigenfunction.
- Utilities
  - [x] Modified Lentz method for continued fractions
  - [ ] WKB approximation for general potentials

Some examples and tutorials are available in the `examples` folder. These are [Pluto.jl](https://plutojl.org/) notebooks.

- [x] Scalar QNMs of RN AdS black hole using the ultraspherical spectral method.

### Time Domain

- **Spherical Symmetry**

  - [ ] Regge-Wheeler-Zerilli equation
    - [ ] Hyperboloidal coordinates
    - [ ] Kerr-Schild coordinates
    - [ ] Tortoise coordinates
  - [ ] Klein-Gordon equation
    - [ ] Hyperboloidal coordinates
    - [ ] Kerr-Schild coordinates
    - [ ] Tortoise coordinates
  - [ ] Einstein equations with a scalar field
    - [ ] Z4 formulation
    - [ ] Hyperboloidal coordinates

- **Axisymmetry**

  - [ ] Teukolsky equation
    - [ ] Hyperboloidal coordinates
    - [ ] Kerr-Schild coordinates
    - [ ] Tortoise coordinates
  - [ ] Klein-Gordon equation
    - [ ] Hyperboloidal coordinates
    - [ ] Kerr-Schild coordinates
    - [ ] Tortoise coordinates

While the `ChebyshevSuite` is inspired by algorithms from [Chebfun](https://www.chebfun.org/), it has been significantly enhanced in this package to improve performance and support arbitrary-precision calculations.

## Acknowledgments

This package is inspired by the following projects, and some of the algorithms are translated / adapted from them:

- [Chebfun](https://www.chebfun.org/)
- [JuliaApproximation/ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl)
- [duetosymmetry/qnm](https://github.com/duetosymmetry/qnm)
- [SciML/MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl)
- [Neal/xsum](https://gitlab.com/radfordneal/xsum)
- [JuliaMath/KahanSummation.jl](https://github.com/JuliaMath/KahanSummation.jl)
- [tomtrogdon/URCMethod.jl](https://github.com/tomtrogdon/URCMethod.jl)

The author would like to thank the developers of these projects for their contributions to the scientific computing community.

## Other Related Projects

- [lucass-carneiro/QuasinormalModes.jl: A Julia package for computing discrete eigenvalues of second-order ODEs](https://github.com/lucass-carneiro/QuasinormalModes.jl): This package focuses on computing quasinormal modes using the `Asymptotic Iteration Method`. However, as far as I know, no one has proven that the `Asymptotic Iteration Method` is guaranteed to converge. So why not use more reliable methods? In any case, it is good to have one more method to compare against.
- [JLRipley314/TeukolskyQNMFunctions.jl: Computes quasinormal modes and eigenfunctions of the Teukolsky equation in HPHC coordinates.](https://github.com/JLRipley314/TeukolskyQNMFunctions.jl): Justin Ripley has developed a package that computes the quasinormal modes of the Teukolsky equation based on his own paper. We also follow his paper for the Teukolsky equation in hyperboloidal coordinates.
