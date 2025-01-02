# PDESuite

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AuroraDysis.github.io/PDESuite.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AuroraDysis.github.io/PDESuite.jl/dev/)
[![Build Status](https://github.com/AuroraDysis/PDESuite.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AuroraDysis/PDESuite.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AuroraDysis/PDESuite.jl/graph/badge.svg?token=C99DVUUULL)](https://codecov.io/gh/AuroraDysis/PDESuite.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Description

This package is a **high-performance** suite designed to compute **arbitrary-precision** solutions for partial differential equations (PDEs) in general relativity. It includes the following numerical methods:

- Chebyshev collocation at Chebyshev points of the first and second kinds
- Finite difference method
- Rectangular collocation
- Ultraspherical spectral method [WIP]
- Ultraspherical rectangular collocation [TODO]

Although the ChebyshevSuite is inspired by algorithms from [Chebfun](https://github.com/chebfun/chebfun), it has been extensively enhanced in this package to optimize performance and enable arbitrary-precision calculations.
