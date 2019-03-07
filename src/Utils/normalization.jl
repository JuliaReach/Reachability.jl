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

"""
    distribute_initial_set(system::InitialValueProblem{<:HybridSystem, <:LazySet)

Distribute the set of initial states to each mode of a hybrid system.

### Input

- `system` -- an initial value problem wrapping a mathematical system (hybrid)
              and a set of initial states

### Output

A new initial value problem with the same hybrid system but where the set of initial
states is the list of tuples `(state, X0)`, for each state in the hybrid system.
"""
function distribute_initial_set(system::InitialValueProblem{<:HybridSystem, <:LazySet})
    HS, X0 = system.s, system.x0
    initial_states = [(loc, X0) for loc in states(HS)]
    return InitialValueProblem(HS, initial_states)
end
