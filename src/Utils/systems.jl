using Systems, LazySets

export ContinuousSystem,
       DiscreteSystem,
       add_dimension,
       next_set

# name alises
const LCS = LinearContinuousSystem
const LDS = LinearDiscreteSystem
const CLCCS = ConstrainedLinearControlContinuousSystem
const CLCDS = ConstrainedLinearControlDiscreteSystem

import Base: *

*(M::AbstractMatrix, input::ConstantInput) =  ConstantInput(M * input.U)

# no input: x' = Ax, x(0) = X0
ContinuousSystem(A::AbstractMatrix, X0::LazySet) = IVP(LCS(A), X0)
DiscreteSystem(A::AbstractMatrix, X0::LazySet) = IVP(LDS(A), X0)

# with constant input: x' = Ax + u, x(0) = X0, u(t) = U
ContinuousSystem(A::AbstractMatrix, X0::LazySet, U::ConstantInput) = IVP(CLCCS(A, eye(A), nothing, U), X0)
ContinuousSystem(A::AbstractMatrix, X0::LazySet, U::LazySet) = ContinuousSystem(A, X0, ConstantInput(U))

DiscreteSystem(A::AbstractMatrix, X0::LazySet, U::ConstantInput) = IVP(CLCDS(A, eye(A), nothing, U), X0)
DiscreteSystem(A::AbstractMatrix, X0::LazySet, U::LazySet) = DiscreteSystem(A, X0, ConstantInput(U))

# with time-varying input: x' = Ax + u, x(0) = X0, u(t) = U(t)
ContinuousSystem(A::AbstractMatrix, X0::LazySet, U::VaryingInput) = IVP(CLCCS(A, eye(A), nothing, U), X0)
ContinuousSystem(A::AbstractMatrix, X0::LazySet, U::Vector{<:LazySet}) = ContinuousSystem(A, X0, VaryingInput(U))

DiscreteSystem(A::AbstractMatrix, X0::LazySet, U::VaryingInput) = IVP(CLCDS(A, eye(A), nothing, U), X0)
DiscreteSystem(A::AbstractMatrix, X0::LazySet, U::Vector{<:LazySet}) = DiscreteSystem(A, X0, VaryingInput(U))

# convenience functions
next_set(inputs::ConstantInput) = next(inputs, 1)[1]
next_set(inputs::AbstractInput, state::Int64) = next(inputs, state)[1]

"""
    add_dimension(A::AbstractMatrix)

Adds an extra zero row and column to a matrix.

### Input

- `A` -- matrix

### Examples

```julia
julia> A = [0.4 0.25; 0.46 -0.67]
2×2 Array{Float64,2}:
 0.4    0.25
 0.46  -0.67
julia> add_dimension(A)
3×3 Array{Float64,2}:
 0.4    0.25  0.0
 0.46  -0.67  0.0
 0.0    0.0   0.0
```
"""
function add_dimension(A::AbstractMatrix)
    n = size(A, 1)
    return vcat(hcat(A, zeros(n)), zeros(n+1).')
end

"""
    add_dimension(X::LazySet)::LazySet

Adds an extra dimension to a LazySet, usually through a Cartesian product.

### Input

- `X` -- a lazy set

### Examples

```jldoctest
julia> X = BallInf(ones(9), 0.5);

julia> dim(X)
9

julia> Xext = add_dimension(X);

julia> dim(Xext)
10

julia> X = ZeroSet(4);

julia> dim(add_dimension(X))
5

julia> typeof(X)
LazySets.ZeroSet{Float64}
```
"""
function add_dimension(X::LazySet)::LazySet
    if X isa ZeroSet
        return ZeroSet(dim(X)+1)
    else
        return X * ZeroSet(1)
    end
end

"""
    add_dimension(cont_sys)

Adds an extra dimension to a continuous system.

### Input

- `cs` -- continuous system

### Examples

```julia
julia> A = sprandn(3, 3, 0.5);
julia> X0 = BallInf(zeros(3), 1.0);
julia> s = ContinuousSystem(A, X0);
julia> sext = add_dimension(s);
julia> statedim(sext)
4
```

If there is an input set, it is also extended:

```julia
julia> U = Ball2(ones(3), 0.1);
julia> s = ContinuousSystem(A, X0, U);
julia> sext = add_dimension(s);
julia> statedim(sext)
4
julia> dim(next_set(inputset(sext)))
4
```

Extending a system with a varying input set:

If there is an input set, it is also extended:

```julia
julia> U = [Ball2(ones(3), 0.1 * i) for i in 1:3];
julia> s = ContinuousSystem(A, X0, U);
julia> sext = add_dimension(s);
julia> statedim(sext)
4
julia> dim(next_set(inputset(sext), 1))
4
```
"""
function add_dimension(cs)
    Aext = add_dimension(cs.s.A)
    X0ext = add_dimension(cs.x0)
    if :U in fieldnames(cs.s)
        Uext = map(add_dimension, inputset(cs))
        return ContinuousSystem(Aext, X0ext, Uext)
    else
        return ContinuousSystem(Aext, X0ext)
    end
end
