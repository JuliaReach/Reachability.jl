# Reachability.jl

[![Build Status](https://travis-ci.org/JuliaReach/Reachability.jl.svg?branch=master)](https://travis-ci.org/JuliaReach/Reachability.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliareach.github.io/Reachability.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/JuliaReach/Reachability.jl/blob/master/LICENSE)
[![Code coverage](http://codecov.io/github/JuliaReach/Reachability.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaReach/Reachability.jl?branch=master)
[![Join the chat at https://gitter.im/JuliaReach/Lobby](https://badges.gitter.im/JuliaReach/Lobby.svg)](https://gitter.im/JuliaReach/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

`Reachability` is a software for reachability analysis and safety property checking that performs flowpipe computation of dynamical systems given by ordinary differential equations models (ODEs) in continuous or discrete time. 
It is written in [Julia](http://julialang.org), a modern high-performance language
for scientific computing.

Currently this package implements algorithms for flowpipe approximation of:

- linear and nonlinear purely continuous or purely discrete ODEs with nondeterministic inputs
- hybrid dynamical systems

## Resources

- [Manual](http://juliareach.github.io/Reachability.jl/dev/)
- [Contributing](https://juliareach.github.io/Reachability.jl/dev/about/#Contributing-1)
- [Release notes of tagged versions](https://github.com/JuliaReach/Reachability.jl/releases)
- [Release log of the development version](https://github.com/JuliaReach/Reachability.jl/wiki/Release-log-tracker)
- [Publications](https://juliareach.github.io/Reachability.jl/dev/publications/)
- [Citations](https://juliareach.github.io/Reachability.jl/dev/citations/)

## Installing

This package requires Julia v1.0 or later. Refer to the [official documentation](https://julialang.org/downloads)
on how to install and run Julia in your system.

### Dependencies

The set computations depend on the core library [`LazySets.jl`](https://github.com/JuliaReach/LazySets.jl), which is also part of the [JuliaReach](https://github.com/JuliaReach/) framework. `LazySets.jl` exploits the principle of lazy (on-demand) evaluation and uses support functions to represent lazy sets. 

The latest stable release of [`LazySets.jl`](https://github.com/JuliaReach/LazySets.jl) is installed automatically when you install `Reachability.jl` (see the installation instructions below). See the [installation section of `LazySets.jl`](https://juliareach.github.io/LazySets.jl/dev/man/getting_started/) for further details.

### Installation

Depending on your needs, choose an appropriate command from the following list
and enter it in Julia's REPL. To activate the `pkg` mode, type `]` (and to leave it, type `<backspace>`).

#### [Install the latest release version](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Adding-registered-packages-1)

```julia
pkg> add https://github.com/JuliaReach/Reachability.jl
```

#### Install the latest development version

```julia
pkg> add https://github.com/JuliaReach/Reachability.jl#master
```

#### [Clone the package for development](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Developing-packages-1)

```julia
pkg> dev https://github.com/JuliaReach/Reachability.jl
```
