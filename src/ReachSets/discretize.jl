const LDS = LinearDiscreteSystem
const CLCDS = ConstrainedLinearControlDiscreteSystem

@inline I(T, n) = Matrix{eltype(A)}(I, n, n)

"""
    discretize(𝑆, δ; [approximation], [exp_method], [sih_method])

Apply an approximation model to `S` obtaining a discrete initial value problem.

### Input

- `𝑆`             -- initial value problem for a continuous affine ODE with
                     non-deterministic inputs
- `δ`             -- step size
- `approximation` -- the method to compute the approximation model for the
                     discretization, choose among:

    - `"forward"`    -- use forward-time interpolation
    - `"backward"`   -- use backward-time interpolation
    - `"firstorder"` -- use first-order approximation of the ODE
    - `"nobloating"` -- do not bloat the initial states

- `exp_method`  -- (optional, default: `"base"`) the method used to take the matrix
                    exponential of the coefficient matrix, choose among:

    - `"base"` -- the scaling and squaring method implemented in Julia base,
                  see `?exp` for details
    - `"pade"` -- use Pade approximant method to compute matrix exponentials of
                  sparse matrices, implemented in `Expokit`
    - `"lazy"` -- compute a wrapper type around the matrix exponential, i.e. using
                  the lazy implementation `SparseMatrixExp` from `LazySets` and
                  the evaluation of the action of the matrix exponential using the
                  `expmv` implementation from `Expokit`

- `sih_method`  -- (optional, default: `"lazy"`) the method used to take the
                    symmetric interval hull operation, choose among:

    - `"concrete"` -- compute the full symmetric interval hull using the function
                      `symmetric_interval_hull` from `LazySets.Approximations`
    - `"lazy"`     -- compute a wrapper set type around symmetric interval hull
                      in a lazy way using `SymmetricIntervalHull`

### Output

The initial value problem of a discrete system.

### Algorithm

Let ``𝑆 : x' = Ax(t) + u(t)``, ``x(0) ∈ \\mathcal{X}_0``, ``u(t) ∈ U`` be the
given continuous affine ODE `𝑆`, where `U` is the set of non-deterministic inputs
and ``\\mathcal{X}_0`` is the set of initial states. Recall that the system
`𝑆` is called homogeneous whenever `U` is the empty set.

Given a step size ``δ``, this function computes a set, `Ω₀`, that guarantees to
contain all the trajectories of ``𝑆`` starting at any ``x(0) ∈ \\mathcal{X}_0``
and for any input function that satisfies ``u(t) ∈ U``, for any ``t ∈ [0, δ]``.

The initial value problem returned by this function consists of the set `Ω₀`
together with the coefficient matrix ``ϕ = e^{Aδ}`` and a transformed
set of inputs if `U` is non-empty.

In the literature, the method to obtain `Ω₀` is called the *approximation model*
and different alternatives have been proposed. See the argument `approximation`
for available options. For the reference to the original papers, see the docstring
of each method.

In the dense-time case, the transformation described is such that the trajectories
of the given continuous system are included in the computed flowpipe of the
discretized system.

In the discrete-time case, there is no bloating of the initial states and the
input is assumed to remain constant between sampled times. Use the option
`approximation="nobloating"` for this setting.

Several methods to compute the matrix exponential are availabe. Use `exp_method`
to select one. For very large systems (~10000×10000), computing the full matrix
exponential is very expensive hence it is preferable to compute the action
of the matrix exponential over vectors when needed. Use the option
`exp_method="lazy"` for this.
"""
function discretize(𝑆::InitialValueProblem{<:AbstractContinuousSystem},
                    δ::Float64;
                    approximation::String="forward",
                    exp_method::String="base",
                    sih_method::String="lazy")

    if approximation in ["forward", "backward"]
        return _discretize_interpolation(𝑆, δ, approximation=approximation,
                    exp_method=exp_method, sih_method=sih_method)
    elseif approximation == "firstorder"
        return _discretize_firstorder(𝑆, δ, exp_method=exp_method)
    elseif approximation == "nobloating"
        return _discretize_nobloating(𝑆, δ, exp_method=exp_method)
    else
        throw(ArgumentError("the approximation model $approximation is unknown"))
    end
end

"""
    exp_Aδ(A::AbstractMatrix, δ::Float64; [exp_method])

Compute the matrix exponential ``e^{Aδ}``.

### Input

- `A`           -- coefficient matrix
- `δ`           -- step size
- `exp_method`  -- (optional, default: `"base"`) the method used to take the matrix
                    exponential of the coefficient matrix, choose among:

    - `"base"` -- the scaling and squaring method implemented in Julia base,
                  see `?exp` for details
    - `"pade"` -- use Pade approximant method to compute matrix exponentials of
                  sparse matrices, as implemented in `Expokit`
    - `"lazy"` -- compute a wrapper type around the matrix exponential, i.e. using
                  the lazy implementation `SparseMatrixExp` from `LazySets` and
                  evaluation of the action of the matrix exponential using the
                  `expmv` implementation in `Expokit`

### Output

A matrix.
"""
function exp_Aδ(A::AbstractMatrix{Float64}, δ::Float64; exp_method="base")
    if exp_method == "base"
        return expmat(Matrix(A*δ))
    elseif exp_method == "lazy"
        return SparseMatrixExp(A*δ)
    elseif exp_method == "pade"
        return padm(A*δ)
    else
       throw(ArgumentError("the exponentiation method $exp_method is unknown"))
    end
end

"""
    ϕ₁(A, δ; [exp_method])

TODO: Add doctring

### Input

- `A`           -- coefficient matrix
- `δ`           -- step size
- `exp_method`  -- (optional, default: `"base"`) the method used to take the matrix
                    exponential of the coefficient matrix, choose among:

    - `"base"` -- the scaling and squaring method implemented in Julia base,
                  see `?exp` for details
    - `"pade"` -- use Pade approximant method to compute matrix exponentials of
                  sparse matrices, as implemented in `Expokit`
    - `"lazy"` -- compute a wrapper type around the matrix exponential, i.e. using
                  the lazy implementation `SparseMatrixExp` from `LazySets` and
                  evaluation of the action of the matrix exponential using the
                  `expmv` implementation in `Expokit`

### Output

A matrix.
"""
function ϕ₁(A, δ; exp_method="base")
    n = size(A, 1)
    if exp_method == "base"
        P = expmat(Matrix([A*δ     sparse(δ*I, n, n)  spzeros(n, n);
                   spzeros(n, 2*n) sparse(δ*I, n, n);
                   spzeros(n, 3*n)]))
        ϕ₁_Aδ = P[1:n, (n+1):2*n]

    elseif exp_method == "lazy"
        P = SparseMatrixExp([A*δ sparse(δ*I, n, n) spzeros(n, n);
                             spzeros(n, 2*n) sparse(δ*I, n, n);
                             spzeros(n, 3*n)])
        ϕ₁_Aδ = sparse(get_columns(P, (n+1):2*n)[1:n, :])

    elseif exp_method == "pade"
        P = padm([A*δ sparse(δ*I, n, n) spzeros(n, n);
                  spzeros(n, 2*n) sparse(δ*I, n, n);
                  spzeros(n, 3*n)])
       ϕ₁_Aδ = P[1:n, (n+1):2*n]

    else
       throw(ArgumentError("the exponentiation method $exp_method is unknown"))
    end

     return ϕ₁_Aδ
end

"""
    ϕ₂(A, δ; [exp_method])

TODO: Add doctring

### Input

- `A`           -- coefficient matrix
- `δ`           -- step size
- `exp_method`  -- (optional, default: `"base"`) the method used to take the matrix
                    exponential of the coefficient matrix, choose among:

    - `"base"` -- the scaling and squaring method implemented in Julia base,
                  see `?exp` for details
    - `"pade"` -- use Pade approximant method to compute matrix exponentials of
                  sparse matrices, as implemented in `Expokit`
    - `"lazy"` -- compute a wrapper type around the matrix exponential, i.e. using
                  the lazy implementation `SparseMatrixExp` from `LazySets` and
                  evaluation of the action of the matrix exponential using the
                  `expmv` implementation in `Expokit`

### Output

A matrix.
"""
function ϕ₂(A, δ; exp_method="base")
    n = size(A, 1)
    if exp_method == "base"
        P = expmat(Matrix([A*δ sparse(δ*I, n, n) spzeros(n, n);
                   spzeros(n, 2*n) sparse(δ*I, n, n);
                   spzeros(n, 3*n)]))
        ϕ₂_Aδ = P[1:n, (2*n+1):3*n]

    elseif exp_method == "lazy"
        P = SparseMatrixExp([A*δ sparse(δ*I, n, n) spzeros(n, n);
                             spzeros(n, 2*n) sparse(δ*I, n, n);
                             spzeros(n, 3*n)])
        ϕ₂_Aδ = sparse(get_columns(P, (2*n+1):3*n)[1:n, :])

    elseif exp_method == "pade"
        P = padm([A*δ sparse(δ*I, n, n) spzeros(n, n);
                  spzeros(n, 2*n) sparse(δ*I, n, n);
                  spzeros(n, 3*n)])
        ϕ₂_Aδ = P[1:n, (2*n+1):3*n]

    else
       throw(ArgumentError("the exponentiation method $exp_method is unknown"))
    end

    return ϕ₂_Aδ
end

"""
    _discretize_firstorder(𝑆, δ; [p], [exp_method])

Apply a first-order approximation model to `S` obtaining a discrete initial value problem.

### Input

- `𝑆`           -- initial value problem for a continuous affine ODE with
                   non-deterministic inputs
- `δ`           -- step size
- `p`           -- (optional, default: `Inf`) parameter in the considered norm
- `exp_method`  -- (optional, default: `base`) the method used to take the matrix
                   exponential of the coefficient matrix, choose among:

    - `base` -- the scaling and squaring method implemented in Julia base,
                see `?exp` for details
    - `pade` -- use Pade approximant method to compute matrix exponentials of
                sparse matrices, as implemented in `Expokit`
    - `lazy` -- compute a wrapper type around the matrix exponential, i.e. using
                the lazy implementation `SparseMatrixExp` from `LazySets` and
                evaluation of the action of the matrix exponential using the
                `expmv` implementation in `Expokit`

### Output

The initial value problem for a discrete system.

### Algorithm

Let us define some notation. Let ``𝑆 : x' = Ax(t) + u(t)``,
``x(0) ∈ \\mathcal{X}_0``, ``u(t) ∈ U`` be the given continuous affine ODE `𝑆`,
where `U` is the set of non-deterministic inputs and ``\\mathcal{X}_0`` is the set
of initial states.

Let ``R_{\\mathcal{X}_0} = \\max_{x ∈ \\mathcal{X}_0} ‖x‖``,
`D_{\\mathcal{X}_0} = \\max_{x, y ∈ \\mathcal{X}_0} ‖x-y‖`` and
``R_{V} = \\max_{u ∈ U} ‖u‖``.

Let ``Ω₀`` be the set defined as:
```math
Ω₀ = ConvexHull(\\mathcal{X}_0, e^{δA}\\mathcal{X}_0 ⊕ δU ⊕ αB_p)
```
where ``α = (e^{δ ‖A‖} - 1 - δ‖A‖)*R_{\\mathcal{X}_0} + R_{U} / ‖A‖)`` and ``B_p`` denotes
the unit ball for the considered norm.

It is proved in [Lemma 1, 1] that the set of states reachable by ``S`` in the time
interval ``[0, δ]``, that we denote ``R_{[0,δ]}(\\mathcal{X}_0)``,
is included in ``Ω₀``:

```math
R_{[0,δ]}(\\mathcal{X}_0) ⊆ Ω₀.
```

Moreover, if `d_H(A, B)` denotes the Hausdorff distance between the sets ``A``
and ``B`` in ``\\mathbb{R}^n``, then

```math
d_H(Ω₀, R_{[0,δ]}(\\mathcal{X}_0)) ≤ \\frac{1}{4}(e^{δ ‖A‖} - 1) D_{\\mathcal{X}_0} + 2α.
```

### Notes

In this implementation, the infinity norm is used by default. To use other norms
substitute `BallInf` with the ball in the appropriate norm. However, note that
not all norms are supported; see the documentation of `?norm` in `LazySets` for
details.

See also [`discr_bloat_interpolation`](@ref) for an alternative algorithm that
uses less conservative bounds.

[1] Le Guernic, C., & Girard, A., 2010, *Reachability analysis of linear systems
using support functions. Nonlinear Analysis: Hybrid Systems, 4(2), 250-262.*
"""
function _discretize_firstorder(𝑆::InitialValueProblem,
                                δ::Float64;
                                p::Float64=Inf,
                                exp_method::String="base")

    # unwrap coefficient matrix and initial states
    A, X0 = 𝑆.s.A, 𝑆.x0 

    # system size; A is assumed square
    n = size(A, 1)

    Anorm = norm(Matrix(A), p)
    RX0 = norm(X0, p)

    # compute exp(A*δ)
    ϕ = exp_Aδ(A, δ, exp_method)

    if islinear(𝑆) # inputdim(𝑆) == 0
        α = (exp(δ*Anorm) - 1. - δ*Anorm) * RX0
        □ = Ballp(p, zeros(n), α)
        Ω0 = ConvexHull(X0, ϕ * X0 ⊕ □)
        return IVP(LDS(ϕ), Ω0)
    elseif isaffine(𝑆)
        Uset = inputset(𝑆)
        if Uset isa ConstantInput
            U = next_set(Uset)
            RU = norm(U, Inf)
            α = (exp(δ*Anorm) - 1.0 - δ*Anorm)*(RX0 + RU/Anorm)
            β = (exp(δ*Anorm) - 1.0 - δ*Anorm)*RU/Anorm
            □α = Ballp(p, zeros(n), α)
            □β = Ballp(p, zeros(n), β)
            Ω0 = ConvexHull(X0, ϕ * X0 ⊕ δ * U + □α)
            Ud = map(u -> δ*u ⊕ □β, U)
            return IVP(CLCDS(ϕ, I(typeof(A), n), nothing, Ud), Ω0)

        elseif Uset isa VaryingInput
            Ud = Vector{LazySet}(undef, length(Uset)) # TODO: concrete type of Uset
            for (i, Ui) in enumerate(Uset)
                RU = norm(Ui, p)
                α = (exp(δ*Anorm) - 1.0 - δ*Anorm)*(RX0 + RU/Anorm)
                β = (exp(δ*Anorm) - 1.0 - δ*Anorm)*RU/Anorm
                □α = Ballp(p, zeros(n), α)
                □β = Ballp(p, zeros(n), β)
                Ω0 = ConvexHull(X0, ϕ * X0 ⊕ δ * Ui ⊕ □α)
                Ud[i] =  δ * Ui ⊕ □β
            end
            Ud = VaryingInput(Ud)
            return IVP(CLCDS(ϕ, I(typeof(ϕ), n), nothing, Ud), Ω0)
        end
    else
        throw(ArgumentError("this function only applies to linear or affine systems"))
    end
end

"""
    _discretize_nobloating(𝑆, δ; [exp_method])

Discretize a continuous system without bloating of the initial states, suitable
for discrete-time reachability.

## Input

- `𝑆`          -- a continuous system
- `δ`          -- step size
- `exp_method` -- (optional, default: `"base"`) the method used to take the matrix
                   exponential of the coefficient matrix, choose among:

    - `"base"` -- the scaling and squaring method implemented in Julia base,
                  see `?exp` for details
    - `"pade"` -- use Pade approximant method to compute matrix exponentials of
                  sparse matrices, as implemented in `Expokit`
    - `"lazy"` -- compute a wrapper type around the matrix exponential, i.e. using
                  the lazy implementation `SparseMatrixExp` from `LazySets` and
                  evaluation of the action of the matrix exponential using the
                  `expmv` implementation in `Expokit`

## Output

A discrete system.

## Algorithm

The transformation implemented here is the following:

- `A -> Phi := exp(A*delta)`
- `U -> V := M*U`
- `X0 -> X0hat := X0`

where `M` corresponds to `Phi1(A, delta)` in Eq. (8) of *SpaceEx: Scalable
Verification of Hybrid Systems.*

In particular, there is no bloating, i.e. we don't bloat the initial states and
dont multiply the input by the step size δ, as required for the dense time case.
"""
function  _discretize_nobloating(𝑆::InitialValueProblem{<:AbstractContinuousSystem},
                                 δ::Float64;
                                 exp_method::String="base")

    # unrwap coefficient matrix and initial states
    A, X0 = 𝑆.s.A, 𝑆.x0

    # compute matrix ϕ = exp(Aδ)
    ϕ = exp_Aδ(A, δ, lazy_expm, pade_expm)

    # early return for homogeneous systems
    if islinear(𝑆)
        Ω0 = X0
        return IVP(LDS(ϕ), Ω0)
    end

    U = inputset(𝑆)
    inputs = next_set(U, 1)

    # compute matrix to transform the inputs
    Phi1Adelta = ϕ₁(A, δ, exp_method)

    discretized_U = Phi1Adelta * inputs

    Ω0 = X0

    if U isa ConstantInput
        return DiscreteSystem(ϕ, Ω0, discretized_U)
    else
        discretized_U = VaryingInput([Phi1Adelta * Ui for Ui in U])
        return DiscreteSystem(ϕ, Ω0, discretized_U)
    end
end

"""
    _discretize_interpolation(𝑆, δ, [approximation], [exp_method], [sih_method])

Compute bloating factors using forward or backward interpolation.

## Input

- `cs`            -- a continuous system
- `δ`             -- step size
- `approximation` -- choose the approximation model among `"forward"` and
                     `"backward"`
- `exp_method`    -- (optional, default: `"base"`) the method used to take the matrix
                     exponential of the coefficient matrix, choose among:

    - `"base"`    -- the scaling and squaring method implemented in Julia base,
                     see `?exp` for details
    - `"pade"`    -- use Pade approximant method to compute matrix exponentials of
                     sparse matrices, as implemented in `Expokit`
    - `"lazy"`    -- compute a wrapper type around the matrix exponential, i.e. using
                     the lazy implementation `SparseMatrixExp` from `LazySets` and
                     evaluation of the action of the matrix exponential using the
                     `expmv` implementation in `Expokit`

- `sih_method`    -- (optional, default: `"lazy"`) the method used to take the
                     symmetric interval hull operation, choose among:

    - `"concrete"` -- compute the full symmetric interval hull
    - `"lazy"`     -- compute a wrapper set type around symmetric interval hull in a
                      lazy way

## Algorithm

See Frehse et al., CAV'11, *SpaceEx: Scalable Verification of Hybrid Systems*,
Lemma 3.

Note that in the unlikely case that A is invertible, the result can also
be obtained directly, as a function of the inverse of A and `e^{At} - I`.

The matrix `P` is such that: `ϕAabs = P[1:n, 1:n]`,
`Phi1Aabsdelta = P[1:n, (n+1):2*n]`, and `Phi2Aabs = P[1:n, (2*n+1):3*n]`.
"""
function _discretize_interpolation(𝑆::InitialValueProblem{<:AbstractContinuousSystem},
                                   δ::Float64;
                                   approximation::String="forward",
                                   exp_method::String="base",
                                   sih_method::String="lazy")

    if sih_method == "lazy"
        sih = SymmetricIntervalHull
    elseif sih_method == "concrete"
        sih = symmetric_interval_hull
    else
        throw(ArgumentError("the method $sih_method is unknown"))
    end

    # unrwap coefficient matrix and initial states
    A, X0 = 𝑆.s.A, 𝑆.x0

    # compute matrix ϕ = exp(Aδ)
    ϕ = exp_Aδ(A, δ, lazy_expm, pade_expm)

    # early return for homogeneous systems
    if islinear(𝑆)
        Ω0 = ConvexHull(X0, ϕ * X0 ⊕ E)
        return IVP(LDS(ϕ), Ω0)
    end
    U = inputset(𝑆)
    inputs = next_set(U, 1)

    # compute the transformation matrix to bloat the initial states
    Phi2Aabs = ϕ₂_Aδ(abs.(A), δ, exp_method=exp_method)

    if isa(inputs, ZeroSet)
        if approximation == "forward" || approximation == "backward"
            Ω0 = ConvexHull(X0, ϕ * X0 + δ * inputs)
        end
    else
        EPsi = sih(Phi2Aabs * sih(A * inputs))
        discretized_U = δ * inputs + EPsi
        if approximation == "forward"
            EOmegaPlus = sih(Phi2Aabs * sih((A * A) * X0))
            Ω0 = ConvexHull(X0, ϕ * X0 + discretized_U + EOmegaPlus)
        elseif approximation == "backward"
            EOmegaMinus = sih(Phi2Aabs * sih((A * A * ϕ) * X0))
            Ω0 = ConvexHull(X0, ϕ * X0 + discretized_U + EOmegaMinus)
        end
    end

    if U isa ConstantInput
        return DiscreteSystem(ϕ, Ω0, discretized_U)
    else
        discretized_U = [δ * Ui + sih(Phi2Aabs * sih(A * Ui)) for Ui in U]
        return DiscreteSystem(ϕ, Ω0, discretized_U)
    end
end
