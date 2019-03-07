"""
    normalize(system)

Normalize a continuous system.

### Input

- `system` -- continuous system

### Output

Either the same system if it already conforms to our required structure, or a
new system otherwise.

### Algorithm

We apply the following normalization steps.

* [`normalize_inputs`](@ref)
"""
function normalize(system::InitialValueProblem{<:Union{AbstractContinuousSystem,
                                                       AbstractDiscreteSystem}})
    return normalize_inputs(system)
end

"""
    normalize_inputs(system)

Normalize the inputs of a continuous system.

### Input

- `system` -- continuous system

### Output

Either the same system if the inputs are of type `AbstractInput`, or a new
system that wraps the inputs in a `ConstantInput`.
"""
function normalize_inputs(system)
    inputdim(system) == 0 && return system
    U = inputset(system)
    U isa AbstractInput && return system
    if !(U isa LazySet)
        throw(ArgumentError("inputs of type $(typeof(U)) are not supported"))
    end
    return IVP(wrap_inputs(system.s, U), system.x0)
end

function wrap_inputs(system::CLCCS, U::LazySet)
    return CLCCS(system.A, system.B, system.X, ConstantInput(U))
end

function wrap_inputs(system::CLCDS, U::LazySet)
    return CLCDS(system.A, system.B, system.X, ConstantInput(U))
end

function wrap_inputs(system::CACCS, U::LazySet)
    return CACCS(system.A, system.B, system.c, system.X, ConstantInput(U))
end

function wrap_inputs(system::CACDS, U::LazySet)
    return CACDS(system.A, system.B, system.c, system.X, ConstantInput(U))
end
