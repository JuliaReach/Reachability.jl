const LDS = LinearDiscreteSystem
const CLCDS = ConstrainedLinearControlDiscreteSystem

@inline Id(n) = Matrix(1.0I, n, n)

"""
    discretize(𝑆, δ; [algorithm], [exp_method], [sih_method], [set_operations])

Apply an approximation model to `S` obtaining a discrete initial value problem.

### Input

- `𝑆`             -- initial value problem for a continuous affine ODE with
                     non-deterministic inputs
- `δ`             -- step size
- `algorithm`     -- (optional, default: `"forward"`) the algorithm used to
                     compute the approximation model for the discretization,
                     choose among:

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

- `sih_method`  -- (optional, default: `"concrete"`) the method used to take the
                    symmetric interval hull operation, choose among:

    - `"concrete"` -- compute the full symmetric interval hull using the function
                      `symmetric_interval_hull` from `LazySets.Approximations`
    - `"lazy"`     -- compute a wrapper set type around symmetric interval hull
                      in a lazy way using `SymmetricIntervalHull`

- `set_operations`  -- (optional, default: `"lazy"`) set operations used for the
                       discretized initial states and transformed inputs, choose among:

    - `"lazy"`     -- use lazy convex hull for the initial states and lazy linear
                      map for the inputs
    - `"zonotope"` -- use concrete zonotope operations (linear map and Minkowski sum)
                      and return zonotopes for both the initial states and the
                      inputs of the discretized system

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
and different alternatives have been proposed. See the argument `algorithm`
for available options. For the reference to the original papers, see the
docstring of each method `discretize_...`.

In the dense-time case, the transformation is such that the trajectories
of the given continuous system are included in the computed flowpipe of the
discretized system.

In the discrete-time case, there is no bloating of the initial states and the
input is assumed to remain constant between sampled times. Use the option
`algorithm="nobloating"` for this setting.

Several methods to compute the matrix exponential are availabe. Use `exp_method`
to select one. For very large systems, computing the full matrix exponential is
expensive hence it is preferable to compute the action of the matrix exponential
over vectors when needed, `e^{δA} v` for each `v`. Use the option
`exp_method="lazy"` for this purpose.
"""
function discretize(𝑆::InitialValueProblem{<:AbstractContinuousSystem},
                    δ::Float64;
                    algorithm::String="forward",
                    exp_method::String="base",
                    sih_method::String="concrete",
                    set_operations::String="lazy")

    if algorithm in ["forward", "backward"]
        return discretize_interpolation(𝑆, δ, algorithm=algorithm,
                    exp_method=exp_method, sih_method=sih_method, set_operations=set_operations)

    elseif algorithm == "firstorder"
        if set_operations != "lazy"
            throw(ArgumentError("the algorithm $algorithm with set operations=$set_operations is not implemented"))
        end
        return discretize_firstorder(𝑆, δ, exp_method=exp_method)

    elseif algorithm == "nobloating"
        if set_operations != "lazy"
            throw(ArgumentError("the algorithm $algorithm with set operations=$set_operations is not implemented"))
        end
        return discretize_nobloating(𝑆, δ, exp_method=exp_method)
    else
        throw(ArgumentError("the algorithm $algorithm is unknown"))
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
                  sparse matrices, implemented in `Expokit`
    - `"lazy"` -- compute a wrapper type around the matrix exponential, i.e. using
                  the lazy implementation `SparseMatrixExp` from `LazySets` and
                  the evaluation of the action of the matrix exponential using the
                  `expmv` implementation from `Expokit`

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

@inline function Pmatrix(A, δ, n)
    return [A*δ     sparse(δ*I, n, n)  spzeros(n, n)    ;
            spzeros(n, 2*n          )  sparse(δ*I, n, n);
            spzeros(n, 3*n          )                   ]
end

"""
    Φ₁(A, δ; [exp_method])

Compute the series

```math
Φ₁(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+1}}{(i+1)!}A^i,
```
where ``A`` is a square matrix of order ``n`` and ``δ ∈ \\mathbb{R}_{≥0}``.

### Input

- `A`           -- coefficient matrix
- `δ`           -- step size
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
### Output

A matrix.

### Algorithm

We use the method from [1]. If ``A`` is invertible, ``Φ₁`` can be computed as

```math
Φ₁(A, δ) = A^{-1}(e^{δA} - I_n).
```

In the general case, implemented in this function, it can be computed as
submatrices of the block matrix

```math
P = \\exp \\begin{pmatrix}
Aδ && δI_n && 0 \\
0 && 0 && δI_n \\
0 && 0 && 0
\\end{array}.
```
It can be shown that `Φ₁(A, δ) = P[1:n, (n+1):2*n]`.

[1] Frehse, Goran, et al. "SpaceEx: Scalable verification of hybrid systems."
International Conference on Computer Aided Verification. Springer, Berlin,
Heidelberg, 2011.
"""
function Φ₁(A, δ; exp_method="base")
    n = size(A, 1)
    if exp_method == "base"
        P = expmat(Matrix(Pmatrix(A, δ, n)))
        Φ₁_Aδ = P[1:n, (n+1):2*n]

    elseif exp_method == "lazy"
        P = SparseMatrixExp(Pmatrix(A, δ, n))
        Φ₁_Aδ = sparse(get_columns(P, (n+1):2*n)[1:n, :])

    elseif exp_method == "pade"
        P = padm(Pmatrix(A, δ, n))
        Φ₁_Aδ = P[1:n, (n+1):2*n]

    else
        throw(ArgumentError("the exponentiation method $exp_method is unknown"))
    end

    return Φ₁_Aδ
end

"""
    Φ₂(A, δ; [exp_method])

Compute the series

```math
Φ₂(A, δ) = ∑_{i=0}^∞ \\dfrac{δ^{i+2}}{(i+2)!}A^i,
```
where ``A`` is a square matrix of order ``n`` and ``δ ∈ \\mathbb{R}_{≥0}``.

### Input

- `A`           -- coefficient matrix
- `δ`           -- step size
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

### Output

A matrix.

### Algorithm

We use the method from [1]. If ``A`` is invertible, ``Φ₂`` can be computed as

```math
Φ₂(A, δ) = A^{-2}(e^{δA} - I_n - δA).
```

In the general case, implemented in this function, it can be computed as
submatrices of the block matrix

```math
P = \\exp \\begin{pmatrix}
Aδ && δI_n && 0 \\
0 && 0 && δI_n \\
0 && 0 && 0
\\end{array}.
```
It can be shown that `Φ₂(A, δ) = P[1:n, (2*n+1):3*n]`.

[1] Frehse, Goran, et al. "SpaceEx: Scalable verification of hybrid systems."
International Conference on Computer Aided Verification. Springer, Berlin,
Heidelberg, 2011.
"""
function Φ₂(A, δ; exp_method="base")
    n = size(A, 1)
    if exp_method == "base"
        P = expmat(Matrix(Pmatrix(A, δ, n)))
        Φ₂_Aδ = P[1:n, (2*n+1):3*n]

    elseif exp_method == "lazy"
        P = SparseMatrixExp(Pmatrix(A, δ, n))
        Φ₂_Aδ = sparse(get_columns(P, (2*n+1):3*n)[1:n, :])

    elseif exp_method == "pade"
        P = padm(Pmatrix(A, δ, n))
        Φ₂_Aδ = P[1:n, (2*n+1):3*n]

    else
       throw(ArgumentError("the exponentiation method $exp_method is unknown"))
    end

    return Φ₂_Aδ
end

"""
    discretize_firstorder(𝑆, δ; [p], [exp_method])

Apply a first-order approximation model to `S` obtaining a discrete initial
value problem.

### Input

- `𝑆`           -- initial value problem for a continuous affine ODE with
                   non-deterministic inputs
- `δ`           -- step size
- `p`           -- (optional, default: `Inf`) parameter in the considered norm
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

### Output

The initial value problem for a discrete system. In particular:

- if the input system is homogeneous, a `LinearDiscreteSystem` is returned,
- otherwise a `ConstrainedLinearControlDiscreteSystem` is returned.

### Algorithm

Let us define some notation. Let

```math
𝑆 : x' = Ax(t) + u(t)
```
with ``x(0) ∈ \\mathcal{X}_0``, ``u(t) ∈ U`` be the given continuous affine ODE
`𝑆`, where `U` is the set of non-deterministic inputs and ``\\mathcal{X}_0``
is the set of initial states.

Define ``R_{\\mathcal{X}_0} = \\max_{x ∈ \\mathcal{X}_0} ‖x‖``,
`D_{\\mathcal{X}_0} = \\max_{x, y ∈ \\mathcal{X}_0} ‖x-y‖`` and
``R_{U} = \\max_{u ∈ U} ‖u‖``. If only the support functions of ``\\mathcal{X}_0``
and ``U`` are known, these values might be hard to compute for some norms. See
`Notes` below for details.

Let ``Ω₀`` be the set defined as:
```math
Ω₀ = ConvexHull(\\mathcal{X}_0, e^{δA}\\mathcal{X}_0 ⊕ δU ⊕ αB_p)
```
where ``α = (e^{δ ‖A‖} - 1 - δ‖A‖)*(R_{\\mathcal{X}_0} + R_{U} / ‖A‖)`` and
``B_p`` denotes the unit ball for the considered ``p``-norm.

It is proved in [Lemma 1, 1] that the set of states reachable by ``S`` in the time
interval ``[0, δ]``, which we denote ``R_{[0,δ]}(\\mathcal{X}_0)`` here,
is included in ``Ω₀``:

```math
R_{[0,δ]}(\\mathcal{X}_0) ⊆ Ω₀.
```

Moreover, if ``d_H(A, B)`` denotes the Hausdorff distance between the sets ``A``
and ``B`` in ``\\mathbb{R}^n``, then

```math
d_H(Ω₀, R_{[0,δ]}(\\mathcal{X}_0)) ≤ \\frac{1}{4}(e^{δ ‖A‖} - 1) D_{\\mathcal{X}_0} + 2α.
```
Hence, the approximation error can be made arbitrarily small by choosing ``δ``
small enough.

Here we allow ``U`` to be a sequence of time varying non-deterministic inputs.

### Notes

In this implementation, the infinity norm is used by default. Other usual norms
are ``p=1`` and ``p=2``. However, note that not all norms are supported; see the
documentation of `?norm` in `LazySets` for the supported norms.

See also [`discretize_interpolation`](@ref) for an alternative algorithm that
uses less conservative bounds.

[1] Le Guernic, C., & Girard, A., 2010, *Reachability analysis of linear systems
using support functions. Nonlinear Analysis: Hybrid Systems, 4(2), 250-262.*
"""
function discretize_firstorder(𝑆::InitialValueProblem,
                               δ::Float64;
                               p::Float64=Inf,
                               exp_method::String="base")

    # unwrap coefficient matrix and initial states
    A, X0 = 𝑆.s.A, 𝑆.x0 

    # system size; A is assumed square
    n = size(A, 1)

    Anorm = norm(Matrix(A), p)
    κ = exp(δ*Anorm) - 1.0 - δ*Anorm
    RX0 = norm(X0, p)

    # compute exp(A*δ)
    ϕ = exp_Aδ(A, δ, exp_method=exp_method)

    if inputdim(𝑆) == 0
        α = κ * RX0
        □α = Ballp(p, zeros(n), α)
        Ω0 = ConvexHull(X0, ϕ * X0 ⊕ □α)
        return IVP(LDS(ϕ), Ω0)

    elseif isaffine(𝑆)
        Uset = inputset(𝑆)
        if Uset isa ConstantInput
            U = next_set(Uset)
            RU = norm(U, p)

            # bloating coefficients
            α = κ*(RX0 + RU/Anorm)
            β = κ*RU/Anorm

            # transformation of the initial states
            □α = Ballp(p, zeros(n), α)
            Ω0 = ConvexHull(X0, ϕ*X0 ⊕ δ*U ⊕ □α)

            # transformation of the inputs
            □β = Ballp(p, zeros(n), β)
            Ud = ConstantInput(δ*U ⊕ □β)
            return IVP(CLCDS(ϕ, Id(n), nothing, Ud), Ω0)

        elseif Uset isa VaryingInput
            Ud = Vector{LazySet}(undef, length(Uset))
            for (i, Ui) in enumerate(Uset)
                RU = norm(Ui, p)

                # bloating factors
                α = κ*(RX0 + RU/Anorm)
                β = κ*RU/Anorm

                if i == 1
                    # transform initial states
                    □α = Ballp(p, zeros(n), α)
                    Ω0 = ConvexHull(X0, ϕ * X0 ⊕ δ * Ui ⊕ □α)
                end

                # transform inputs
                □β = Ballp(p, zeros(n), β)
                Ud[i] = δ*Ui ⊕ □β
            end
            Ud = VaryingInput(Ud)
            return IVP(CLCDS(ϕ, Id(n), nothing, Ud), Ω0)
        else
            throw(ArgumentError("input of type $(typeof(U)) is not allowed"))
        end
    else
        throw(ArgumentError("this function only applies to linear or affine systems"))
    end
end

"""
    discretize_nobloating(𝑆, δ; [exp_method])

Discretize a continuous system without bloating of the initial states, suitable
for discrete-time reachability.

### Input

- `𝑆`          -- a continuous system
- `δ`          -- step size
- `exp_method` -- (optional, default: `"base"`) the method used to take the matrix
                  exponential of the coefficient matrix, choose among:

    - `"base"` -- the scaling and squaring method implemented in Julia base,
                  see `?exp` for details
    - `"pade"` -- use Pade approximant method to compute matrix exponentials of
                  sparse matrices, implemented in `Expokit`
    - `"lazy"` -- compute a wrapper type around the matrix exponential, i.e. using
                  the lazy implementation `SparseMatrixExp` from `LazySets` and
                  the evaluation of the action of the matrix exponential using the
                  `expmv` implementation from `Expokit`

### Output

The initial value problem for a discrete system. In particular:

- if the input system is homogeneous, a `LinearDiscreteSystem` is returned,
- otherwise a `ConstrainedLinearControlDiscreteSystem` is returned.

### Algorithm

Let us define some notation. Let

```math
𝑆 : x' = Ax(t) + u(t)
```
with ``x(0) ∈ \\mathcal{X}_0``, ``u(t) ∈ U`` be the given continuous affine ODE
`𝑆`, where `U` is the set of non-deterministic inputs and ``\\mathcal{X}_0``
is the set of initial states.

The approximation model implemented in this function is such that there is no bloating,
i.e. we don't bloat the initial states and don't multiply the input by the step
size δ, as required for the dense time case.

The transformations are:

- ``Φ ← \\exp^{Aδ}``
- ``Ω₀ ← \\mathcal{X}_0``
- ``V ← Φ₁(A, δ)U(k)``, where ``Φ₁(A, δ)`` is defined in
  [`Φ₁(A, δ; [exp_method])`](@ref).

Here we allow ``U`` to be a sequence of time varying non-deterministic input sets.
"""
function  discretize_nobloating(𝑆::InitialValueProblem{<:AbstractContinuousSystem},
                                δ::Float64;
                                exp_method::String="base")

    # unwrap coefficient matrix and initial states
    A, X0 = 𝑆.s.A, 𝑆.x0

    # compute matrix ϕ = exp(Aδ)
    ϕ = exp_Aδ(A, δ, exp_method=exp_method)

    # initial states remain unchanged
    Ω0 = copy(X0)

    # early return for homogeneous systems
    if inputdim(𝑆) == 0
        return IVP(LDS(ϕ), Ω0)
    end

    # compute matrix to transform the inputs
    Phi1Adelta = Φ₁(A, δ, exp_method=exp_method)
    U = inputset(𝑆)
    Ud = map(ui -> Phi1Adelta * ui, U)

    return IVP(CLCDS(ϕ, Id(size(A, 1)), nothing, Ud), Ω0)
end

"""
    discretize_interpolation(𝑆, δ; [algorithm], [exp_method], [sih_method])

Compute bloating factors using forward or backward interpolation.

### Input

- `𝑆`             -- a continuous system
- `δ`             -- step size
- `algorithm`     -- choose the algorithm to compute the approximation model
                     among `"forward"` and `"backward"`
- `exp_method`    -- (optional, default: `"base"`) the method used to take the matrix
                     exponential of the coefficient matrix, choose among:

    - `"base"` -- the scaling and squaring method implemented in Julia base,
                  see `?exp` for details
    - `"pade"` -- use Pade approximant method to compute matrix exponentials of
                  sparse matrices, implemented in `Expokit`
    - `"lazy"` -- compute a wrapper type around the matrix exponential, i.e. using
                  the lazy implementation `SparseMatrixExp` from `LazySets` and
                  the evaluation of the action of the matrix exponential using the
                  `expmv` implementation from `Expokit`

- `sih_method`    -- (optional, default: `"concrete"`) the method used to take the
                     symmetric interval hull operation, choose among:

    - `"concrete"` -- compute the full symmetric interval hull
    - `"lazy"`     -- compute a wrapper set type around symmetric interval hull in a
                      lazy way

### Output

The initial value problem for a discrete system. In particular:

- if the input system is homogeneous, a `LinearDiscreteSystem` is returned,
- otherwise a `ConstrainedLinearControlDiscreteSystem` is returned.

## Algorithm

Let us define some notation. Let

```math
𝑆 : x' = Ax(t) + u(t)
```
with ``x(0) ∈ \\mathcal{X}_0``, ``u(t) ∈ U`` be the given continuous affine ODE
`𝑆`, where `U` is the set of non-deterministic inputs and ``\\mathcal{X}_0``
is the set of initial states.

The transformations are:

- ``Φ ← \\exp^{Aδ}``,
- ``Ω₀ ← ConvexHull(\\mathcal{X}_0, Φ\\mathcal{X}_0 ⊕ δU(0) ⊕ Eψ(U(0), δ) ⊕ E^+(\\mathcal{X}_0, δ))``,
- ``V ← δU(k) ⊕ Eψ(U(k), δ)``.

Here we allow ``U`` to be a sequence of time varying non-deterministic input sets.

For the definition of the sets ``Eψ`` and ``E^+`` see [1]. The  "backward" method
uses ``E^-``.

[1] Frehse, Goran, et al. "SpaceEx: Scalable verification of hybrid systems."
International Conference on Computer Aided Verification. Springer, Berlin,
Heidelberg, 2011.
"""
function discretize_interpolation(𝑆::InitialValueProblem{<:AbstractContinuousSystem},
                                   δ::Float64;
                                   algorithm::String="forward",
                                   exp_method::String="base",
                                   sih_method::String="concrete",
                                   set_operations::String="lazy")

    # used to dispatch on the value of the set operation
    set_operations = Symbol(set_operations)

    if sih_method == "concrete"
        sih = symmetric_interval_hull
    elseif sih_method == "lazy"
        sih = SymmetricIntervalHull
    else
        throw(ArgumentError("the method $sih_method is unknown"))
    end

    # unwrap coefficient matrix and initial states
    A, X0 = 𝑆.s.A, 𝑆.x0

    # compute matrix ϕ = exp(Aδ)
    ϕ = exp_Aδ(A, δ, exp_method=exp_method)

    # compute the transformation matrix to bloat the initial states
    Phi2Aabs = Φ₂(abs.(A), δ, exp_method=exp_method)

    if algorithm == "forward"
        Einit = sih(Phi2Aabs * sih((A * A) * X0))     # use E⁺
    elseif algorithm == "backward"
        Einit = sih(Phi2Aabs * sih((A * A * ϕ) * X0)) # use E⁻
    else
        throw(ArgumentError("the algorithm $approximation is unknown"))
    end

    # early return for homogeneous systems
    if inputdim(𝑆) == 0
        Ω0 = _discretize_interpolation_homog(X0, ϕ, Einit, Val(set_operations))
        return IVP(LDS(ϕ), Ω0)
    end

    U = inputset(𝑆)
    U0 = next_set(U, 1)

    Eψ0 = sih(Phi2Aabs * sih(A * U0))
    Ω0, Ud = _discretize_interpolation_inhomog(δ, U0, U, X0, ϕ, Einit, Eψ0, Phi2Aabs, Val(set_operations))

    return IVP(CLCDS(ϕ, Id(size(A, 1)), nothing, Ud), Ω0)
end

# version using lazy sets and operations
function _discretize_interpolation_homog(X0, ϕ, Einit, set_operations::Val{:lazy})
    Ω0 = ConvexHull(X0, ϕ * X0 ⊕ Einit)
    return Ω0
end

# version using lazy sets and operations
function _discretize_interpolation_inhomog(δ, U0, U, X0, ϕ, Einit, Eψ0, Phi2Aabs, set_operations::Val{:lazy})
    Ω0 = ConvexHull(X0, ϕ * X0 ⊕ δ*U0 ⊕ Eψ0 ⊕ Einit)

    if U isa ConstantInput
        Ud = ConstantInput(δ*U0 ⊕ Eψ0)

    elseif U isa VaryingInput
        Ud = Vector{LazySet}(undef, length(U))
        for (k, Uk) in enumerate(U)
            Eψk = sih(Phi2Aabs * sih(A * Uk))
            Ud[k] = δ * Uk ⊕ Eψk
        end
        Ud = VaryingInput(Ud)

    else
        throw(ArgumentError("input of type $(typeof(U)) is not allowed"))
    end
    return Ω0, Ud
end

# version using concrete operations with zonotopes
function _discretize_interpolation_homog(X0, ϕ, Einit, set_operations::Val{:zonotope})
    Einit = convert(Zonotope, Einit)
    Z1 = convert(Zonotope, X0)
    Z2 = linear_map(ϕ, minkowski_sum(X0, Einit))
    Ω0 = overapproximate(ConvexHull(Z1, Z2), Zonotope)
    return Ω0
end

# version using concrete operations with zonotopes
function _discretize_interpolation_inhomog(δ, U0, U, X0, ϕ, Einit, Eψ0, Phi2Aabs, set_operations::Val{:zonotope})
    Einit = convert(Zonotope, Einit)
    Eψ0 = convert(Zonotope, Eψ0)
    Z1 = X0
    ϕX0 = linear_map(ϕ, X0)
    U0 = convert(Zonotope, U0)
    δI = Matrix(δ*I, size(U0))
    δU0 = linear_map(δI, U0)
    Z2 = reduce(minkowski_sum, [ϕX0, δU0, Eψ0, Einit])
    Ω0 = overapproximate(ConvexHull(Z1, Z2), Zonotope)

    if U isa ConstantInput
        Ud = ConstantInput(minkowski_sum(δU0, Eψ0))

    elseif U isa VaryingInput
        Ud = Vector{LazySet}(undef, length(U))
        for (k, Uk) in enumerate(U)
            Eψk = convert(Zonotope, sih(Phi2Aabs * sih(A * Uk)))
            Uk = convert(Zonotope, Uk)
            δUk = linear_map(δI, Uk)
            Ud[k] = minkowski_sum(δUk, Eψk)
        end
        Ud = VaryingInput(Ud)

    else
        throw(ArgumentError("input of type $(typeof(U)) is not allwed"))
    end
    return Ω0, Ud
end
