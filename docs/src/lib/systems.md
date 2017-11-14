# Systems

This module provides types to represent systems of affine ODEs in both discrete
and continuous time.

```@contents
Pages = ["systems.md"]
Depth = 3
```

```@meta
CurrentModule = Reachability.Systems
```

## Different types of systems

Every system inherits from `AbstractSystem`.

```@docs
AbstractSystem
```

We support the following two concrete types of systems. 

### Discrete system

A discrete system consists of a matrix representing the system dynamics, a set
of initial states, a set of nondeterministic inputs, and a discretization step
δ.

```@docs
DiscreteSystem
dim(S::DiscreteSystem)
```

### Continuous system

A continuous system consists of a matrix representing the system dynamics, a set
of initial states, and a set of nondeterministic inputs.

```@docs
ContinuousSystem
dim(S::ContinuousSystem)
```

## Nondeterministic inputs

The above systems may contain nondeterministic inputs, which are wrapped in
special types. Every nondeterministic input representation inherits from
`NonDeterministicInput`.

```@docs
NonDeterministicInput
```

The inputs are closely related to a [`DiscreteSystem`](@ref) in the sense that
for each discrete time step the input set may change. The `InputState` type
allows an iteration through the inputs over time.

```@docs
InputState
```

### Constant nondeterministic inputs

Constant nondeterministic inputs are chosen from a set of values that does not
change over time. Note that, while the set is constant, the inputs themselves
vary over time.

```@docs
ConstantNonDeterministicInput
```

### Time-varying nondeterministic inputs

Time-varying nondeterministic inputs are chosen from a set of values that
changes over time (with each time step).

```@docs
TimeVaryingNonDeterministicInput
```
