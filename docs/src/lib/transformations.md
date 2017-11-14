# Transformations

This module provides functions to apply coordinate transformations to
[Systems](@ref).

```@contents
Pages = ["transformations.md"]
Depth = 3
```

```@meta
CurrentModule = Reachability.Transformations
```

## Interface

This module exports a single function that works as an interface. It dispatches
which transformation to apply using a string argument.

```@docs
transform
```

## Types of transformations

We currently only support a Schur transformation.

```@docs
transform_schur
```
