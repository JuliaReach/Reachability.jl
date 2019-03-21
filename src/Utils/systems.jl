export ContinuousSystem,
       DiscreteSystem,
       add_dimension,
       next_set

# name alises
const LCS = LinearContinuousSystem
const LDS = LinearDiscreteSystem
const CLCS = ConstrainedLinearContinuousSystem
const CLDS = ConstrainedLinearDiscreteSystem
const CLCCS = ConstrainedLinearControlContinuousSystem
const CLCDS = ConstrainedLinearControlDiscreteSystem
const CACCS = ConstrainedAffineControlContinuousSystem
const CACDS = ConstrainedAffineControlDiscreteSystem
const CACS = ConstrainedAffineContinuousSystem
const CADS = ConstrainedAffineDiscreteSystem
const BBCS = BlackBoxContinuousSystem
const CBBCS = ConstrainedBlackBoxContinuousSystem
const CBBCCS = ConstrainedBlackBoxContinuousSystem

export LCS, LDS, CLCS, CLDS, CLCCS, CLCDS, CACCS, CACDS, IVP, BBCS, CBBCS, CBBCCS

import Base: *
import LazySets.constrained_dimensions

*(M::AbstractMatrix, input::ConstantInput) =  ConstantInput(M * input.U)

# convenience functions
next_set(inputs::ConstantInput) = collect(nextinput(inputs, 1))[1]
next_set(inputs::AbstractInput, state::Int64) = collect(nextinput(inputs, state))[1]

"""
    add_dimension(A::AbstractMatrix, m=1)

Adds an extra zero row and column to a matrix.

### Input

- `A` -- matrix
- `m` -- (optional, default: `1`) the number of extra dimensions

### Examples

```jldoctest add_dimension_test
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

```jldoctest add_dimension_test
julia> add_dimension(A, 2)
4×4 Array{Float64,2}:
 0.4    0.25  0.0  0.0
 0.46  -0.67  0.0  0.0
 0.0    0.0   0.0  0.0
 0.0    0.0   0.0  0.0
```
"""
function add_dimension(A::AbstractMatrix, m=1)
    n = size(A, 1)
    return vcat(hcat(A, zeros(n, m)), zeros(m, n+m))
end

"""
    add_dimension(X::LazySet, m=1)::LazySet

Adds an extra dimension to a LazySet through a Cartesian product.

### Input

- `X` -- a lazy set
- `m` -- (optional, default: `1`) the number of extra dimensions

### Examples

```jldoctest add_dimension_set
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

More than one dimension can be added passing the second argument:

```jldoctest add_dimension_set
julia> Xext = add_dimension(BallInf(zeros(10), 0.1), 4);

julia> dim(Xext)
14
```

### Notes

In the special case that the given set is a zero set, instead of cartesian product
a new zero set with extended dimensions is returned.
"""
function add_dimension(X::LazySet, m=1)::LazySet
    return X * ZeroSet(m)
end

function add_dimension(X::ZeroSet, m=1)::ZeroSet
    return ZeroSet(dim(X)+m)
end

"""
    add_dimension(cs, m=1)

Adds an extra dimension to a continuous system.

### Input

- `cs` -- continuous system
- `m` -- (optional, default: `1`) the number of extra dimensions

### Examples

```jldoctest add_dimension_cont_sys
julia> using MathematicalSystems
julia> A = sprandn(3, 3, 0.5);
julia> X0 = BallInf(zeros(3), 1.0);
julia> s = InitialValueProblem(LinearContinuousSystem(A), X0);
julia> sext = add_dimension(s);
julia> statedim(sext)
4
```

If there is an input set, it is also extended:

```jldoctest add_dimension_cont_sys
julia> U = ConstantInput(Ball2(ones(3), 0.1));
julia> s = InitialValueProblem(ConstrainedLinearControlContinuousSystem(A, Matrix(1.0I, size(A)), nothing, U), X0);
julia> sext = add_dimension(s);
julia> statedim(sext)
4
julia> dim(next_set(inputset(sext)))
4
```

Extending a system with a varying input set:

If there is an input set, it is also extended:

```jldoctest add_dimension_cont_sys
julia> U = VaryingInput([Ball2(ones(3), 0.1 * i) for i in 1:3]);
julia> s = InitialValueProblem(ConstrainedLinearControlContinuousSystem(A, Matrix(1.0I, size(A)), nothing, U), X0);
julia> sext = add_dimension(s);
julia> statedim(sext)
4
julia> dim(next_set(inputset(sext), 1))
4
```

Extending a varing input set with more than one extra dimension:

```jldoctest add_dimension_cont_sys
julia> sext = add_dimension(s, 7);
julia> statedim(sext)
10
julia> dim(next_set(inputset(sext), 1))
10
```
"""
function add_dimension(cs, m=1)
    Aext = add_dimension(cs.s.A, m)
    X0ext = add_dimension(cs.x0, m)
    if hasmethod(inputset, Tuple{typeof(cs.s)})
        Uext = map(x -> add_dimension(x, m), inputset(cs))
        s = ConstrainedLinearControlContinuousSystem(Aext, Matrix(1.0I, size(Aext)), nothing, Uext)
    else
        s = LinearContinuousSystem(Aext)
    end
    return InitialValueProblem(s, X0ext)
end

"""
    constrained_dimensions(HS::HybridSystem)::Dict{Int,Vector{Int}}

For each location, compute all dimensions that are constrained in the invariant
or the guard of any outgoing transition.

### Input

- `HS`  -- hybrid system

### Output

A dictionary mapping the index of each location ``ℓ`` to the dimension indices
that are constrained in ``ℓ``.
"""
function constrained_dimensions(HS::HybridSystem)::Dict{Int,Vector{Int}}
    result = Dict{Int,Vector{Int}}()
    sizehint!(result, nstates(HS))
    for mode in states(HS)
        vars = Vector{Int}()
        append!(vars, constrained_dimensions(stateset(HS, mode)))
        for transition in out_transitions(HS, mode)
            append!(vars, constrained_dimensions(stateset(HS, transition)))
        end
        result[mode] = unique(vars)
    end

    return result
end
