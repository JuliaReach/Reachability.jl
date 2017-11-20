# Transformations

This module provides functions to apply coordinate transformations to
[Systems](@ref) using matrix decompositions.

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

## Schur transform

The real Schur decomposition is of the form

$U^TAU = \begin{pmatrix}
T_{11} & T_{12} &\cdots & T_{1b} \\
0 & T_{22} & \cdots & T_{2b} \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & T_{bb}
\end{pmatrix}$
where $T_{ij}$ are 2x2 matrices


```@docs
schur_transform
```

## Examples
