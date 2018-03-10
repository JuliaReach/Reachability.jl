using Systems

export ContinuousSystem, DiscreteSystem

const Atype = Union{AbstractMatrix, SparseMatrixExp}

import Base: *

*(M::AbstractMatrix, input::ConstantInput) =  ConstantInput(M * input.U)

# no input
ContinuousSystem(A::Atype, X0::LazySet) = IVP(LinearContinuousSystem(A), ZeroSet(size(A, 1)))

# with constant input
ContinuousSystem(A::Atype, X0::LazySet, U::LazySet) = IVP(ConstrainedLinearControlContinuousSystem(A, nothing, nothing, ConstantInput(U)), X0)

# with time-varying input
ContinuousSystem(A::Atype, X0::LazySet, U::Vector{<:LazySet}) = IVP(ConstrainedLinearControlContinuousSystem(A, nothing, nothing, VaryingInput(U)), X0)

# no input
DiscreteSystem(A::Atype, X0::LazySet) = IVP(LinearDiscreteSystem(A), ZeroSet(size(A, 1)))

# constructor that creates a ConstantNonDeterministicInput
DiscreteSystem(A::Atype, X0::LazySet, U::LazySet) = IVP(ConstrainedLinearControlDiscreteSystem(A, nothing, nothing, ConstantInput(U)), X0)

# constructor that creates a TimeVaryingNonDeterministicInput
DiscreteSystem(A::Atype, X0::LazySet, U::Vector{<:LazySet}) =IVP(ConstrainedLinearControlDiscreteSystem(A, nothing, nothing, VaryingInput(U)), X0)
