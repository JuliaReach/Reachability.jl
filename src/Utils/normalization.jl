# accepted types of non-deterministic inputs (non-canonical form)
const UNCF = Union{<:LazySet{N}, Vector{<:LazySet{N}}, <:AbstractInput} where {N}

# accepted types for the state constraints (non-canonical form)
const XNCF = Union{<:LazySet{N}, Nothing} where {N}

"""
    normalize(system::AbstractSystem)

Transform a mathematical system to a normalized (or canonical) form.

### Input

- `system` -- system; it can be discrete or continuous

### Output

Either the same system if it already conforms to a canonical form, or a new
system otherwise.

### Notes

The normalization procedure consists of transforming a given system type into a
"canonical" format that is used internally. More details are given below.

### Algorithm

The implementation of `normalize` exploits `MathematicalSystems`'s' types, which
carry information about the problem as a type parameter.

Homogeneous ODEs of the form ``x' = Ax, x ∈ \\mathcal{X}`` are canonical if the
associated problem is a `ConstrainedLinearContinuousSystem` and `A` is a matrix.
This type does not handle non-deterministic inputs.

Note that a `LinearContinuousSystem` does not consider constraints on the
state-space (such as an invariant); to specify state constraints, use a
`ConstrainedLinearContinuousSystem`. If the passed system is a `LinearContinuousSystem`
(i.e. no constraints) then the normalization fixes a universal set (`Universe`)
as the constraint set.

The generalization to canonical systems with constraints and possibly time-varying
non-deterministic inputs is considered next. These systems are of the form
``x' = Ax + u, u ∈ \\mathcal{U}, x ∈ \\mathcal{X}``. The system type is
`ConstrainedLinearControlContinuousSystem`, where `A` is a matrix, `X` is a set
and `U` is an input, that is, any concrete subtype of `AbstractInput`.

If `U` is not given as an input, normalization accepts either a `LazySet`, or
a vector of `LazySet`s. In these cases, the sets are wrapped around an appropriate
concrete input type.

If the system does not conform to a canonical form, the implementation tries
to make the transformation; otherwise an error is thrown. In particular, ODEs
of the form ``x' = Ax + Bu`` are mapped into ``x' = Ax + u, u ∈ B\\mathcal{U}``,
where now ``u`` has the same dimensions as ``x``.

The transformations described above are analogous in the discrete case, i.e.
``x_{k+1} = A x_k`` and ``x_{k+1} = Ax_{k} + u_k, u_k ∈ \\mathcal{U}, x_k ∈ \\mathcal{X}``
for the linear and affine cases respectively.
"""
function normalize(system::AbstractSystem)
    throw(ArgumentError("the system type $(typeof(system)) is currently not supported"))
end

# x' = Ax, in the continuous case
# x+ = Ax, in the discrete case
function normalize(system::LCS{N, AN}) where {N, AN<:AbstractMatrix{N}}
    n = statedim(system)
    X = Universe(n)
    return CLCS(system.A, X)
end

function normalize(system::LDS{N, AN}) where {N, AN<:AbstractMatrix{N}}
    n = statedim(system)
    X = Universe(n)
    return CLDS(system.A, X)
end

# x' = Ax, x ∈ X in the continuous case
# x+ = Ax, x ∈ X in the discrete case
for CL_S in (:CLCS, :CLDS)
    @eval begin
        function normalize(system::$CL_S{N, AN, ST}) where {N, AN<:AbstractMatrix{N}, ST<:LazySet{N}}
            n = statedim(system)
            X = _wrap_invariant(stateset(system), n)
            return $CL_S(system.A, X)
        end
    end
end

# x' = Ax + u, x ∈ X, u ∈ U in the continuous case
# x+ = Ax + u, x ∈ X, u ∈ U in the discrete case
for CLC_X in (:CLCCS, :CLCDS)
    @eval begin
        function normalize(system::$CLC_X{N, IdentityMultiple{N}, ST}) where {N, ST<:LazySet{N}}
            λ = system.B.M.λ
            if λ == oneunit(λ)
                return system
            else
                n = statedim(system)
                X = _wrap_invariant(stateset(system), n)
                U = _wrap_inputs(system.U, Diagonal(fill(λ, n)))
                return $CLC_X(system.A, I(n, N), X, U)
            end
        end
    end
end

# x' = Ax + Bu, x ∈ X, u ∈ U in the continuous case
# x+ = Ax + Bu, x ∈ X, u ∈ U in the discrete case
for CLC_S in (:CLCCS, :CLCDS)
    @eval begin
        function normalize(system::$CLC_S{N, AN, BN, XT, UT}) where {N, AN<:AbstractMatrix{N}, BN<:AbstractMatrix{N}, XT<:XNCF, UT<:UNCF}
            n = statedim(system)
            X = _wrap_invariant(stateset(system), n)
            U = _wrap_inputs(system.U, system.B)
            $CLC_S(system.A, I(n, N), X, U)
        end
    end
end

# x' = Ax + Bu + c, x ∈ X, u ∈ U in the continuous case
# x+ = Ax + Bu + c, x ∈ X, u ∈ U in the discrete case
for CAC_S in (:CACCS, :CACDS)
    @eval begin
        function normalize(system::$CAC_S{N, AN, BN, VN, XT, UT}) where {N, AN<:AbstractMatrix{N}, BN<:AbstractMatrix{N}, VN<:AbstractVector{N}, XT<:XNCF, UT<:UNCF}
            n = statedim(system)
            X = _wrap_invariant(system.X, n)
            U = _wrap_inputs(system.U, system.B, system.c)
            $CAC_S(system.A, I(n, N), X, U)
        end
    end
end

_wrap_invariant(X::LazySet, n::Int) = X
_wrap_invariant(X::Nothing, n::Int) = Universe(n)

_wrap_inputs(U::AbstractInput, B::IdentityMultiple) = U
_wrap_inputs(U::AbstractInput, B::AbstractMatrix) = map(u -> B*u, U)
_wrap_inputs(U::LazySet, B::IdentityMultiple) = ConstantInput(U)
_wrap_inputs(U::LazySet, B::AbstractMatrix) = ConstantInput(B*U)
_wrap_inputs(U::Vector{<:LazySet}, B::IdentityMultiple) = VaryingInput(U)
_wrap_inputs(U::Vector{<:LazySet}, B::AbstractMatrix) = VaryingInput(map(u -> B*u, U))

_wrap_inputs(U::AbstractInput, B::IdentityMultiple, c::AbstractVector) = U ⊕ c
_wrap_inputs(U::AbstractInput, B::AbstractMatrix, c::AbstractVector) = map(u -> B*u ⊕ c, U)
_wrap_inputs(U::LazySet, B::IdentityMultiple, c::AbstractVector) = ConstantInput(U ⊕ c)
_wrap_inputs(U::LazySet, B::AbstractMatrix, c::AbstractVector) = ConstantInput(B*U ⊕ c)
_wrap_inputs(U::Vector{<:LazySet}, B::IdentityMultiple, c::AbstractVector) = VaryingInput(map(u -> u ⊕ c, U))
_wrap_inputs(U::Vector{<:LazySet}, B::AbstractMatrix, c::AbstractVector) = VaryingInput(map(u -> B*u ⊕ c, U))

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
