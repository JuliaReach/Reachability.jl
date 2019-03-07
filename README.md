# Reachability.jl

[![Build Status](https://travis-ci.org/JuliaReach/Reachability.jl.svg?branch=master)](https://travis-ci.org/JuliaReach/Reachability.jl)
[![Docs latest](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliareach.github.io/Reachability.jl/latest/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/JuliaReach/Reachability.jl/blob/master/LICENSE)
[![Code coverage](http://codecov.io/github/JuliaReach/Reachability.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaReach/Reachability.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

`Reachability` is a software for reachability analysis and safety property checking that performs flowpipe computation of dynamical systems given by ordinary differential equations models (ODEs) in continuous or discrete time. 
It is written in [Julia](http://julialang.org), a modern high-performance language
for scientific computing.

Currently this package implements algorithms that can handle:

- flowpipe computation of affine ODEs with nondeterministic inputs
- hybrid dynamical systems (hybrid automata) with nondeterministic affine ODEs in each mode

## Resources

- [Manual](http://juliareach.github.io/Reachability.jl/latest/)
- [Contributing](https://juliareach.github.io/Reachability.jl/latest/about.html#Contributing-1)
- [Release notes of tagged versions](https://github.com/JuliaReach/Reachability.jl/releases)
- [Release log of the development version](https://github.com/JuliaReach/Reachability.jl/wiki/Release-log-tracker)
- [Publications](https://juliareach.github.io/Reachability.jl/latest/publications.html)


## Installing

This package requires Julia v1.0 or later. Refer to the [official documentation](https://julialang.org/downloads)
on how to install and run Julia in your system.

### Dependencies

The set computations depend on the core library [LazySets.jl](https://github.com/JuliaReach/LazySets.jl), which is also part of the [JuliaReach](https://github.com/JuliaReach/) framework. `LazySets` exploits the principle of lazy (on-demand) evaluation and uses support functions to represent lazy sets. 

The latest stable release of [LazySets.jl](https://github.com/JuliaReach/LazySets.jl) is installed automatically when you install `Reachability.jl`, see installation instructions below. You can always install the development version via `Pkg.clone`; see the [installation section of LazySets.jl](https://juliareach.github.io/LazySets.jl/latest/man/getting_started.html#Setup-1) for further details.

### Installation

To install this package, use the following command inside Julia's REPL:
```julia
using Pkg
Pkg.clone("https://github.com/JuliaReach/Reachability.jl")
```

## Updating

To checkout the latest version, do
```julia
Pkg.checkout("Reachability")
````
