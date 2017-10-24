# Reachability.jl

[![Build Status](https://travis-ci.org/JuliaReach/Reachability.jl.svg?branch=master)](https://travis-ci.org/JuliaReach/Reachability.jl)
[![Docs latest](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliareach.github.io/Reachability.jl/latest/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/JuliaReach/Reachability.jl/blob/master/LICENSE.md)

`Reachability` is a software for reach set approximation and safety properties
of affine ordinary differential equations (ODEs) with nondeterminsitic inputs.

It is written in [Julia](http://julialang.org), a modern high-performance language
for scientific computing.

## Resources

- [Manual](http://juliareach.github.io/Reachability.jl/latest/)
- [Contributing](http://juliareach.github.io/Reachability.jl/latest/about.html)

## Installing

### Dependencies

This package requires Julia v0.6 or later. Refer to the official documentation on how to install it for your system. 

To install the [LazySets](https://github.com/JuliaReach/LazySets.jl) dependency, follow the instructions in the [installation section](https://github.com/JuliaReach/LazySets.jl#installing).

### Installation

To install this package, use the following command inside Julia's REPL:
```julia
Pkg.clone("https://github.com/JuliaReach/Reachability.jl")
```

## Updating

To checkout the latest version, do 
```julia
Pkg.checkout("Reachability")
````
