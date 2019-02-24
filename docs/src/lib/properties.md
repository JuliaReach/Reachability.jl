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
check(::Conjunction, ::LazySet)
Disjunction
check(::Disjunction, ::LazySet)
```

### Specific properties


```@docs
SafeStatesProperty
check(::SafeStatesProperty, ::LazySet)
BadStatesProperty
check(::BadStatesProperty, ::LazySet)
```
