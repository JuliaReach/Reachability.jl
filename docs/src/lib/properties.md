# Properties

This module provides representations of (safety) properties.

```@contents
Pages = ["properties.md"]
Depth = 3
```

```@meta
CurrentModule = Reachability.Properties
```

## General property interface

```@docs
Property
```

### Boolean combination of properties

```@docs
Conjunction
dim(::Conjunction)
check(::Conjunction, ::LazySet{N}) where {N<:Real}
Disjunction
dim(::Disjunction)
check(::Disjunction, ::LazySet{N}) where {N<:Real}
```

### Specific properties


```@docs
SafeStatesProperty
dim(::SafeStatesProperty)
check(::SafeStatesProperty, ::LazySet)
BadStatesProperty
dim(::BadStatesProperty)
check(::BadStatesProperty, ::LazySet)
```
