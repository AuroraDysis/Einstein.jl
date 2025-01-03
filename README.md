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

This package is a **high-performance** suite designed to compute **arbitrary-precision** solutions for partial differential equations (PDEs) in general relativity. It includes the following numerical methods:

- Chebyshev collocation at Chebyshev points of the first and second kinds
- Finite difference method
- Rectangular collocation
- Ultraspherical spectral method
- Ultraspherical rectangular collocation [TODO]

Although the ChebSuite is inspired by algorithms from [Chebfun](https://github.com/chebfun/chebfun), it has been extensively enhanced in this package to optimize performance and enable arbitrary-precision calculations.
