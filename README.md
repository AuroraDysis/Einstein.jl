# PDESuite

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AuroraDysis.github.io/PDESuite.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AuroraDysis.github.io/PDESuite.jl/dev/)
[![Build Status](https://github.com/AuroraDysis/PDESuite.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AuroraDysis/PDESuite.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AuroraDysis/PDESuite.jl/graph/badge.svg?token=C99DVUUULL)](https://codecov.io/gh/AuroraDysis/PDESuite.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

> [!WARNING]
> The package is currently under development and is not yet ready for production use. I am migrating the code from a private repository to this public one. The documentation is also being updated to reflect the changes made to the codebase. Please be patient and check back later for updates.

## Description

This package is a high-performance suite designed to compute arbitrary-precision solutions for partial differential equations (PDEs) in general relativity. It includes the following numerical methods:

- Chebyshev collocation at Chebyshev points of the first and second kinds [WIP]  
- Finite difference method [WIP]  
- Rectangular collocation at Chebyshev points of the second kind (TODO)  
- Ultraspherical method (TODO)  
- Ultraspherical rectangular pseudospectral method (TODO)  

Although the ChebyshevSuite is inspired by algorithms from [Chebfun](https://github.com/chebfun/chebfun), it has been extensively enhanced in this package to optimize performance and enable arbitrary-precision calculations.
