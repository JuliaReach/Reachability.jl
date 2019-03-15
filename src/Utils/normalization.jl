# linear ivp in canonical form x' = Ax
const LCF = LCS{N, AN} where {N, AN<:AbstractMatrix{N}}

# affine ivp with constraints in canonical form x' = Ax + u, u ∈ U, x ∈ X
const ACF = CLCCS{N, AN, IdentityMultiple{N}, XT, UT} where {N, AN<:AbstractMatrix{N}, XT<:LazySet{N}, UT<:AbstractInput}

# type union of canonical forms
const CF = Union{<:LCF, <:ACF}

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
"canonical" format that is used internally. The implementation of `normalize`
exploits `MathematicalSystems`'s' types, which carry information about the problem
as a type parameter.

The type union `CF` defines the initial value problems considered canonical,
i.e. which do not require normalization. More details are given below.

Homogeneous ODEs of the form ``x' = Ax`` are canonical if the associated
problem is a `LinearContinuousSystem` and `A` is a matrix.
Note that it does not consider constraints on the state-space (such as an
invariant); in this case you may have to use a
`ConstrainedLinearControlContinuousSystem`. Moreover, this type does not handle
non-deterministic inputs.

The generalization to canonical systems with constraints and possibly time-varying
non-deterministic inputs is considered next. These systems are of the form
``x' = Ax + u, u ∈ \\mathcal{U}, x ∈ \\mathcal{X}``. The system type is
`ConstrainedLinearControlContinuousSystem`, where `A` is a matrix, `X` is a set
and `U` is an input, that is, any concrete subtype of `AbstractInput`.

If `U` is not given as an input, normalization accepts either as `LazySet`, or
a vector of `LazySet`'s. In these cases, the sets are wrapped around an appropriate
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

# systems already in canonical form, i.e. which don't need normalization
normalize(system::CF) = system

# x' = Ax + Bu, x ∈ X, u ∈ U in the continuous case
# x+ = Ax + Bu, x ∈ X, u ∈ U in the discrete case
for CLCS in (:CLCCS, :CLCDS)
    @eval begin
        function normalize(system::$CLCS{N, AN, BN, XT, UT}) where {N, AN<:AbstractMatrix{N}, BN<:AbstractMatrix{N}, XT<:XNCF, UT<:UNCF}
            n = statedim(system)
            X = _wrap_invariant(system.X, n)
            U = _wrap_inputs(system.U, system.B)
            $CLCS(system.A, I(n, N), X, U)
        end
    end
end

# x' = Ax + Bu + c, x ∈ X, u ∈ U in the continuous case
# x+ = Ax + Bu + c, x ∈ X, u ∈ U in the discrete case
for CACS in (:CACCS, :CACDS)
    @eval begin
        function normalize(system::$CACS{N, AN, BN, VN, XT, UT}) where {N, AN<:AbstractMatrix{N}, BN<:AbstractMatrix{N}, VN<:AbstractVector{N}, XT<:XNCF, UT<:UNCF}
            n = statedim(system)
            X = _wrap_invariant(system.X, n)
            U = _wrap_inputs(system.U, system.B, system.c)
            $CACS(system.A, I(n, N), X, U)
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
