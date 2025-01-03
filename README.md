# GRSuite

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AuroraDysis.github.io/GRSuite.jl/stable/) -->

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AuroraDysis.github.io/GRSuite.jl/dev/)
[![Build Status](https://github.com/AuroraDysis/GRSuite.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AuroraDysis/GRSuite.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AuroraDysis/GRSuite.jl/graph/badge.svg?token=C99DVUUULL)](https://codecov.io/gh/AuroraDysis/GRSuite.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

> [!WARNING]  
> I'm still actively developing it and migrating code from private repositories to this one, so it may currently lack some features.

## Motivation and Purpose

After developing hundreds of PDE solvers for various projects, I decided to create `GRSuite.jl` because I couldn't find a library that was both highly convenient and effective for solving partial differential equations (PDEs) in the context of general relativity (GR).

When solving GR PDEs in spherical symmetry or even axisymmetry, the requirements are often far more demanding than in 3+1 cases. These include the need for very long-time simulations, conservation of quantities, highly accurate late-time behavior, and precise computation of quasinormal modes. While I have addressed these challenges individually in the past, I grew tired of copying and pasting solutions from old repositories. Instead, I decided to develop a comprehensive, feature-rich, and highly efficient library. This is my answer to those challenges.

## Description

**GRSuite.jl** is a **high-performance** suite designed to compute **arbitrary-precision solutions** of PDEs in GR. It is planned to support solving the following built-in equations (some of which were developed for previous projects and will be ported to this library):

- **Spherical Symmetry**

  - [ ] Regge-Wheeler-Zerilli equation in hyperboloidal, Kerr-Schild, and tortoise coordinates
  - [ ] Klein-Gordon equation in hyperboloidal, Kerr-Schild, and tortoise coordinates
  - [ ] Einstein equations with a scalar field using the Z4 formulation

- **Axisymmetry**

  - [ ] Teukolsky equation in hyperboloidal, Kerr-Schild, and tortoise coordinates
  - [ ] Klein-Gordon equation in hyperboloidal, Kerr-Schild, and tortoise coordinates

- **Quasinormal Modes**
  - [ ] Compute quasinormal modes for the Schwarzschild black hole using the Regge-Wheeler-Zerilli equation. Use the continued fraction method to determine the eigenvalue and the ultraspherical spectral method in hyperboloidal coordinates to determine the eigenfunction.
  - [ ] Compute quasinormal modes for the Teukolsky equation using the continued fraction method. Apply the Cook-Zalutskiy spectral approach for the angular sector to determine the eigenvalue and the ultraspherical spectral method in hyperboloidal coordinates to determine the eigenfunction.

`GRSuite.jl` currently includes utilities for the following numerical methods:

- **Spatial Discretization**
  - [x] Chebyshev collocation at Chebyshev points of the first and second kinds (Most of algorithms are translated from Chebfun)
  - [x] Finite difference method
  - [x] Rectangular collocation method (Most of algorithms are translated from Chebfun)
  - [x] Ultraspherical spectral method (Most of algorithms are translated from Chebfun)
  - [ ] Ultraspherical rectangular collocation [TODO]
- **Quasinormal Mode Computation**
  - [x] Modified Lentz method for continued fractions
  - [ ] Cook-Zalutskiy spectral method [TODO]

While the `ChebyshevSuite` is inspired by algorithms from [Chebfun](https://www.chebfun.org/), it has been significantly enhanced in this package to improve performance and support arbitrary-precision calculations.

This project is under active development, so feel free to give it a try!
