__precompile__()
"""
Module to handle systems of affine ODEs with nondeterministic inputs.
"""
module Systems

using LazySets

export AbstractSystem,
       ContinuousSystem,
       DiscreteSystem,
       NonDeterministicInput,
       ConstantNonDeterministicInput,
       TimeVaryingNonDeterministicInput

import Base: *


#=
Nondeterministic inputs
=#


"""
Abstract type representing a nondeterministic input. The input can be either constant
or time-varying. In both cases it is represented by an iterator.
"""
abstract type NonDeterministicInput end


"""
Type that represent the state of a `NonDeterministicInput`.

### Fields
- `set`   -- current set
- `index` -- index in the iteration
"""
struct InputState
    # set representation
    set::LazySet
    # iteration index
    index::Int64

    # default constructor
    InputState(set::LazySet, index::Int64) = new(set, index)
end
InputState(set::LazySet) = InputState(set, 1)


"""
Type that represents a constant nondeterministic input.

The iteration over this set is such that its `state` is a tuple
(`set`, `index`), where `set` is the value of the input, represented as a
`LazySet`, and `index` counts the number of times this iterator was called. Its
length is infinite, since the input is defined for all times. The index of the
input state is always contantly 1.

### Fields

- `U` -- `LazySet`
"""
struct ConstantNonDeterministicInput <: NonDeterministicInput
    # input
    U::LazySet

    ConstantNonDeterministicInput(U::LazySet) = new(U)
end
ConstantNonDeterministicInput() = ConstantNonDeterministicInput(VoidSet(0))


Base.start(NDInput::ConstantNonDeterministicInput) = InputState(NDInput.U)
Base.next(NDInput::ConstantNonDeterministicInput, state) = InputState(NDInput.U)
Base.done(NDInput::ConstantNonDeterministicInput, state) = false
Base.eltype(::Type{ConstantNonDeterministicInput}) = Type{InputState}
Base.length(NDInput::ConstantNonDeterministicInput) = 1


function *(M::AbstractMatrix{Float64}, NDInput::ConstantNonDeterministicInput)
    return ConstantNonDeterministicInput(M * start(NDInput).set)
end


"""
Type that represents a time-varying nondeterministic input.

The iteration over this set is such that its `state` is a tuple
(`set`, `index`), where `set` is the value of the input, represented as an array
of `LazySet`s, and `index` counts the number of times this iterator was called.
Its length corresponds to the number of elements in the given array. The index
of the input state increases from 1 and corresponds at each time to the array
index in the input array.

### Fields

- `U` -- array containing `LazySet`s
"""
struct TimeVaryingNonDeterministicInput <: NonDeterministicInput
    # input sequence
    U::Vector{<:LazySet}

    TimeVaryingNonDeterministicInput(U::Vector{<:LazySet}) = new(U)
end


Base.start(NDInput::TimeVaryingNonDeterministicInput) = InputState(NDInput.U[1])
Base.next(NDInput::TimeVaryingNonDeterministicInput, state) = InputState(NDInput.U[state.index+1], state.index+1)
Base.done(NDInput::TimeVaryingNonDeterministicInput, state) = state.index > length(NDInput.U)
Base.eltype(::Type{TimeVaryingNonDeterministicInput}) = Type{InputState}
Base.length(NDInput::TimeVaryingNonDeterministicInput) = length(NDInput.U)


#=
Systems
=#


"""
Abstract type representing a system of affine ODEs.
"""
abstract type AbstractSystem end


"""
Type that represents a continous-time affine system with nondeterministic inputs,

``x'(t) = Ax(t) + u(t)``,

where:

- ``A`` is a square matrix
- ``x(0) ∈ X0``
- ``u(t) ∈ U(t)``, where ``U(t)`` is a piecewise-constant set-valued function
  defined over ``[t1, t1+δ], ... , [tN, tN+δ]`` for some δ

### Fields

- `A`  -- square matrix
- `X0` -- set of initial states
- `U`  -- nondeterministic inputs
"""
struct ContinuousSystem <: AbstractSystem
    # system's matrix
    A::AbstractMatrix{Float64}
    # initial states
    X0::LazySet
    # nondeterministic inputs
    U::NonDeterministicInput

    # default constructor
    ContinuousSystem(A::AbstractMatrix{Float64},
                     X0::LazySet,
                     U::NonDeterministicInput) =
        new(A, X0, U)
end
ContinuousSystem(A::AbstractMatrix{Float64},
                 X0::LazySet) =
    ContinuousSystem(A, X0, ConstantNonDeterministicInput(VoidSet(size(A, 1))))
ContinuousSystem(A::AbstractMatrix{Float64},
                 X0::LazySet,
                 U::LazySet) =
    ContinuousSystem(A, X0, ConstantNonDeterministicInput(U))
ContinuousSystem(A::AbstractMatrix{Float64},
                 X0::LazySet,
                 U::Array{<:LazySet, 1}) =
    ContinuousSystem(A, X0, TimeVaryingNonDeterministicInput(U))


"""
    dim(S)

Dimension of a continuous system.

### Input

- `S` -- continuous system
"""
function dim(S::ContinuousSystem)
    return size(S.A, 1)
end


"""
Type that represents a discrete-time affine system with nondeterministic inputs,

``x_{k+1} = A x_{k} + u_{k}``,

where

- ``A ``is a square matrix
- ``x(0) ∈ X0``
- ``u_{k} ∈ U_{k}``, where ``U_{k}`` is a piecewise-constant set-valued function
  defined over ``[t1, t1+δ], ... , [tN, tN+δ]`` for some δ

### Fields

- `A`  -- square matrix
- `X0` -- set of initial states
- `U`  -- nondeterministic inputs
- `δ`  -- discretization step
"""
struct DiscreteSystem <: AbstractSystem
    # system's matrix
    A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}}
    # initial states
    X0::LazySet
    # nondeterministic inputs
    U::NonDeterministicInput
    # discretization step
    δ::Float64

    # default constructor with some domain checks
    DiscreteSystem(A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}},
                   X0::LazySet,
                   δ::Float64,
                   U::NonDeterministicInput) =
        (δ < 0.
            ? throw(DomainError())
            : new(A, X0, U, δ))
end
DiscreteSystem(A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}},
               X0::LazySet,
               δ::Float64) =
    DiscreteSystem(A, X0, δ, ConstantNonDeterministicInput(VoidSet(size(A, 1))))
DiscreteSystem(A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}},
               X0::LazySet,
               δ::Float64,
               U::LazySet) =
    DiscreteSystem(A, X0, δ, ConstantNonDeterministicInput(U))
DiscreteSystem(A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}},
               X0::LazySet,
               δ::Float64,
               U::Array{<:LazySet, 1}) =
    DiscreteSystem(A, X0, δ, TimeVaryingNonDeterministicInput(U))


"""
    dim(S)

Dimension of a discrete system.

### Input

- `S` -- discrete system
"""
function dim(S::DiscreteSystem)
    return size(S.A, 1)
end

end  # module
