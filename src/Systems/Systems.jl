__precompile__()
"""
Module to handle systems of affine ODEs with non-deterministic inputs.
"""
module Systems

using LazySets

export System,
       ContinuousSystem,
       DiscreteSystem,
       NonDeterministicInput,
       ConstantNonDeterministicInput,
       TimeVaryingNonDeterministicInput

import Base: *

"""
Abstract type representing a nondeterministic input. The input can be either constant
or time-varying. In both cases it is represented by an iterator.
"""
abstract type NonDeterministicInput end

"""
Type that represent the state of a constant or time-varying input, as a pair of
sf (set represented by its support function), and an index.
"""
struct InputState
    sf::LazySet
    index::Int64

    InputState(sf, index) = new(sf, index)
    InputState(sf) = new(sf, 1)
end

"""
Type that represents a constant non-deterministic input.

The iteration over this set is such that its `state` is a tuple (sf, index), where
sf is the value of the input, represented as a LazySet set, and `index` counts
the number of times this iterator was called. Its length is infinite, since the input
is defined for all times. The index of the input state is always contantly 1.
"""
struct ConstantNonDeterministicInput <: NonDeterministicInput
    # input
    U::LazySet

    ConstantNonDeterministicInput(U) = new(U)
end
ConstantNonDeterministicInput() = ConstantNonDeterministicInput(VoidSet(0))

Base.start(NDInput::ConstantNonDeterministicInput) = InputState(NDInput.U)
Base.next(NDInput::ConstantNonDeterministicInput, state) = InputState(NDInput.U)
Base.done(NDInput::ConstantNonDeterministicInput, state) = false
Base.eltype(::Type{ConstantNonDeterministicInput}) = Type{InputState}
Base.length(NDInput::ConstantNonDeterministicInput) = 1

import Base: *

function *(M::AbstractMatrix{Float64}, NDInput::ConstantNonDeterministicInput)
    return ConstantNonDeterministicInput(M * start(NDInput).sf)
end

"""
Type that represents a time-varying non-deterministic input.

The iteration over this set is such that its `state` is a tuple (sf, index), where
sf is the value of the input, represented as a LazySet set, and `index` counts
the number of times this iterator was called. Its length corresponds to the number
of elements in the given array of support functions. The index of the input state
increases from 1 and corresponds at each time to the array index in the input array U.
"""
struct TimeVaryingNonDeterministicInput <: NonDeterministicInput
    # input sequence
    U::Array{<:LazySet, 1}

    TimeVaryingNonDeterministicInput(U) = new(U)
end

Base.start(NDInput::TimeVaryingNonDeterministicInput) = InputState(NDInput.U[1])
Base.next(NDInput::TimeVaryingNonDeterministicInput, state) = InputState(NDInput.U[state.index+1], state.index+1)
Base.done(NDInput::TimeVaryingNonDeterministicInput, state) = state.index > length(NDInput.U)
Base.eltype(::Type{TimeVaryingNonDeterministicInput}) = Type{InputState}
Base.length(NDInput::TimeVaryingNonDeterministicInput) = length(NDInput.U)

"""
Abstract type representing a system of affine ODEs.
"""
abstract type System end

"""
Type that represents an continous-time affine system with non-deterministic inputs,

x'(t) = Ax(t) + u(t),

where:

- A is a square matrix
- x(0) ∈ X0
- u(t) ∈ U(t), where U(t) is a piecewise-constant set-valued function defined
  over [t1, t1+δ], ... , [tN, tN+δ]
"""
struct ContinuousSystem <: System

    # system's matrix
    A::AbstractMatrix{Float64}

    # initial set of states
    X0::LazySet

    # input
    U::NonDeterministicInput

    ContinuousSystem(A::AbstractMatrix{Float64}, X0::LazySet, U::NonDeterministicInput) = new(A, X0, U)
    ContinuousSystem(A::AbstractMatrix{Float64}, X0::LazySet) = new(A, X0, ConstantNonDeterministicInput(VoidSet(size(A, 1))))
    ContinuousSystem(A::AbstractMatrix{Float64}, X0::LazySet, U::LazySet) = new(A, X0, ConstantNonDeterministicInput(U))
    ContinuousSystem(A::AbstractMatrix{Float64}, X0::LazySet, U::Array{<:LazySet, 1}) = new(A, X0, TimeVaryingNonDeterministicInput(U))
end

# dimension of a System (number of independent variables)
function dim(S::ContinuousSystem)
    return size(S.A, 1)
end

"""
Type that represents a discrete-time affine system with non-deterministic inputs,

x_{k+1} = A x_{k} + u_{k},

where

- A is a square matrix
- x(0) ∈ X0
- u_{k} ∈ U_{k}, where U_{k} is a piecewise-constant set-valued function defined over [t1, t1+δ], ... , [tN, tN+δ]
"""
struct DiscreteSystem <: System
    # system's matrix
    A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}}

    # initial set of states
    X0::LazySet

    # discretization
    δ::Float64

    # input
    U::NonDeterministicInput

    DiscreteSystem(A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}}, X0::LazySet, δ::Float64) = new(A, X0, δ, ConstantNonDeterministicInput(VoidSet(size(A, 1))))
    DiscreteSystem(A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}}, X0::LazySet, δ::Float64, U::LazySet) = new(A, X0, δ, ConstantNonDeterministicInput(U))
    DiscreteSystem(A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}}, X0::LazySet, δ::Float64, U::Array{<:LazySet, 1}) = new(A, X0, δ, TimeVaryingNonDeterministicInput(U))
    DiscreteSystem(A::Union{AbstractMatrix{Float64}, SparseMatrixExp{Float64}}, X0::LazySet, δ::Float64, U::NonDeterministicInput) = new(A, X0, δ, U)
end

# dimension of a System (number of independent variables)
function dim(S::DiscreteSystem)
    return size(S.A, 1)
end

end  # module
