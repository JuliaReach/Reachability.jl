"""
    discretize(𝑆, δ; [approx_model], [expm_method], [sih_method])

Apply an approximation model to `S` obtaining a discrete initial value problem.

### Input

- `𝑆`            -- initial value problem for a continuous affine ODE with
                    nondeterministic inputs
- `δ`            -- step size
- `approx_model` -- the method to compute the approximation model for the
                    discretization, choose among:

    - `forward`    -- use forward-time interpolation
    - `backward`   -- use backward-time interpolation
    - `firstorder` -- use first order approximation of the ODE
    - `nobloating` -- do not bloat the initial states

- `expm_method`  -- (optional, default: `base`) the method used to take the matrix
                    exponential of the coefficients matrix, choose among:

    - `base` -- the scaling and squaring method implemented in Julia base,
                see `?exp` for details
    - `pade` -- use Pade approximant method to compute matrix exponentials of
                sparse matrices, as implemented in `Expokit`
    - `lazy` -- compute a wrapper type around the matrix exponential, i.e. using
                the lazy implementation `SparseMatrixExp` from `LazySets` and
                evaluation of the action of the matrix exponential using the
                `expmv` implementation in `Expokit`

- `sih_method`  -- (optional, default: `lazy`) the method used to take the
                    symmetric interval hull operation, choose among:

    - `concrete` -- compute the full symmetric interval hull
    - `lazy`     -- compute a wrapper set type around symmetric interval hull in a
                    lazy way

### Output

The initial value problem of a discrete system.

### Algorithm

Let ``𝑆 : x' = Ax(t) + u(t)``, ``x(0) ∈ \\mathcal{X}_0``, ``u(t) ∈ U`` be the
given continuous affine ODE `𝑆`, where `U` is the set of nondeterministic inputs
and ``\\mathcal{X}_0`` is the set of initial states. Recall that the system
`𝑆` is called homogeneous whenever `U` is the empty set.

Given a step size ``δ``, this function computes a set, `Ω0`, that guarantees to
contain all the trajectories of ``𝑆`` starting at any ``x(0) ∈ \\mathcal{X}_0``
and for any input function that satisfies ``u(t) ∈ U``, for any ``t ∈ [0, δ]``.

The initial value problem returned by this function consists of the set `Ω0`
together with the coefficients matrix ``ϕ = e^{Aδ}`` and a transformed
set of inputs if `U` is non-empty.

In the literature, the method to obtain `Ω0` is called the *approximation model*
and different alternatives have been proposed. See the argument `approx_model`
for available options. For the reference to the original papers, see the docstring
of each method.

The transformation described above allows to do dense time reachability, i.e.
such that the trajectories of the given continuous system are included in the computed
flowpipe of the discretized system. In particular, if you don't want to consider
bloating, i.e. discrete-time reachability, use `approx_model="nobloating"`.

Several methods to compute the matrix exponential are availabe. Use `expm_method`
to select one. For very large systems (~10000×10000), computing the full matrix
exponential is very expensive hence it is preferable to compute the action
of the matrix exponential over vectors when needed. Use the option
`expm_method="lazy"` for this.
"""
function discretize(𝑆::InitialValueProblem{<:AbstractContinuousSystem},
                    δ::Float64;
                    approx_model::String="forward",
                    expm_method::String="base",
                    sih_method::String="lazy")

    if !isaffine(S)
        throw(ArgumentError("`discretize` is only implemented for affine ODEs"))
    end

    if approx_model in ["forward", "backward"]
        return discr_bloat_interpolation(𝑆, δ, approx_model, expm_method, sih_method)
    elseif approx_model == "firstorder"
        return _discretize_first_order(𝑆, δ, expm_method)
    elseif approx_model == "nobloating"
        return discr_no_bloat(𝑆, δ, expm_method)
    else
        throw(ArgumentError("the approximation model $approx_model is unknown"))
    end
end

function exp_Aδ(A::AbstractMatrix, δ::Float64, expm_method="base")
    if expm_method == "base"
        return expmat(Matrix(A*δ))
    elseif expm_method == "lazy"
        return SparseMatrixExp(A*δ)
    elseif expm_method == "pade"
        return padm(A*δ)
    else
       throw(ArgumentError("the exponentiation method $expm_method is unknown"))
    end
end

"""
    _discretize_first_order(𝑆, δ, [p], [expm_method])

Apply a first order approximation model to `S` obtaining a discrete initial value problem.

### Input

- `𝑆` -- initial value problem for a continuous affine ODE with nondeterministic inputs
- `δ` -- step size
- `p` -- (optional, default: `Inf`) parameter in the considered norm

### Output

The initial value problem for a discrete system.

### Algorithm

Let us define some notation. Let ``𝑆 : x' = Ax(t) + u(t)``,
``x(0) ∈ \\mathcal{X}_0``, ``u(t) ∈ U`` be the given continuous affine ODE `𝑆`,
where `U` is the set of nondeterministic inputs and ``\\mathcal{X}_0`` is the set
of initial states.

Let ``R_{X0} = \\max_{x∈ X0} ||x||`` and `D_{X0} = \\max_{x, y ∈ X0} ||x-y||``
and ``R_{V} = \\max_{u ∈ U} ||u||``.

It is proved [Lemma 1, 1] that the set `Ω0` defined as

```math
Ω0 = CH(\\mathcal{X}_0, e^{δA}\\mathcal{X}_0 ⊕ δU ⊕ αB_p)
```
where ``α = (e^{δ||A||} - 1 - δ||A||)*R_{X0} + R_{U} / ||A||)`` and ``B_p`` denotes
the unit ball for the considered norm, is such that the set of states reachable
by ``S`` in the time interval ``[0, δ]``, that we denote ``R_{[0,δ]}(X0)``,
is included in `Ω0`. Moreover, if `d_H(A, B)` denotes the Hausdorff distance
between the sets ``A`` and ``B`` in ``\\mathbb{R}^n``, then

```math
d_H(Ω0, R_{[0,δ]}(X0)) ≤ \\frac{1}{4}(e^{δ||A||} - 1) D_{X0} + 2α.
```

### Notes

In this implementation, the infinity norm is used by default. To use other norms
substitute `BallInf` with the ball in the appropriate norm. However, note that
not all norms are supported; see the documentation of `?norm` in `LazySets` for
details.

[1] Le Guernic, C., & Girard, A., 2010, *Reachability analysis of linear systems
using support functions. Nonlinear Analysis: Hybrid Systems, 4(2), 250-262.*

### Notes

See also [`discr_bloat_interpolation`](@ref) for an alternative algorithm that
uses less conservative bounds.
"""
function _discretize_first_order(𝑆::InitialValueProblem, δ::Float64,
                                 p::Float64=Inf,
                                 expm_method::String="base")

    # unwrap coeffs matrix and initial states
    A, X0 = 𝑆.s.A, 𝑃.x0 

    # system size; A is assumed square
    n = size(A, 1)

    Anorm = norm(Matrix(A), Inf)
    RX0 = norm(X0, Inf)

    # compute exp(A*δ)
    ϕ = exp_Aδ(A, δ, expm_method)

    if islinear(𝑆) # inputdim(𝑆) == 0
        α = (exp(δ*Anorm) - 1. - δ*Anorm) * RX0
        □ = Ballp(p, zeros(n), α)
        Ω0 = CH(X0, ϕ * X0 + □)
        return InitialValueProblem(LinearDiscreteSystem(ϕ), Ω0))  # @system x' = ϕ*x, x(0) ∈ Ω0
    # ---- TODO ---- seguir aca ----
    else
        # affine case; TODO: unify Constant and Varying input branches?
        Uset = inputset(cont_sys)
        if Uset isa ConstantInput
            U = next_set(Uset)
            RU = norm(U, Inf)
            α = (exp(δ*Anorm) - 1. - δ*Anorm)*(RX0 + RU/Anorm)
            β = (exp(δ*Anorm) - 1. - δ*Anorm)*RU/Anorm
            Ω0 = CH(X0, ϕ * X0 + δ * U + Ball2(zeros(size(ϕ, 1)), α))
            discr_U =  δ * U + Ball2(zeros(size(ϕ, 1)), β)
            return DiscreteSystem(ϕ, Ω0, discr_U)
        elseif Uset isa VaryingInput
            discr_U = Vector{LazySet}(undef, length(Uset))
            for (i, Ui) in enumerate(Uset)
                RU = norm(Ui, Inf)
                α = (exp(δ*Anorm) - 1. - δ*Anorm)*(RX0 + RU/Anorm)
                β = (exp(δ*Anorm) - 1. - δ*Anorm)*RU/Anorm
                Ω0 = CH(X0, ϕ * X0 + δ * Ui + Ball2(zeros(size(ϕ, 1)), α))
                discr_U[i] =  δ * Ui + Ball2(zeros(size(ϕ, 1)), β)
            end
            return DiscreteSystem(ϕ, Ω0, discr_U)
        end
    end
end

"""
    discr_no_bloat(cont_sys, δ, pade_expm, lazy_expm)

Discretize a continuous system without bloating of the initial states, suitable
for discrete-time reachability.

## Input

- `cont_sys`     -- a continuous system
- `δ`            -- step size
- `pade_expm`    -- if `true`, use Pade approximant method to compute the
                    matrix exponential
- `lazy_expm`    -- if `true`, compute the matrix exponential in a lazy way
                    (suitable for very large systems)

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
function discr_no_bloat(cont_sys::InitialValueProblem{<:AbstractContinuousSystem},
                        δ::Float64,
                        pade_expm::Bool,
                        lazy_expm::Bool)

    # unrwap coefficients matrix and initial states
    A, X0 = cont_sys.s.A, cont_sys.x0

    # compute matrix ϕ = exp(Aδ)
    ϕ = exp_Aδ(A, δ, lazy_expm, pade_expm)

    # early return for homogeneous systems
    if cont_sys isa IVP{<:LinearContinuousSystem}
        Ω0 = X0
        return DiscreteSystem(ϕ, Ω0)
    end

    U = inputset(cont_sys)
    inputs = next_set(U, 1)

    n = size(A, 1)

    # compute matrix to transform the inputs
    if lazy_expm
        P = SparseMatrixExp([A*δ sparse(δ*I, n, n) spzeros(n, n);
                             spzeros(n, 2*n) sparse(δ*I, n, n);
                             spzeros(n, 3*n)])
        Phi1Adelta = sparse(get_columns(P, (n+1):2*n)[1:n, :])
    else
        if pade_expm
            P = padm([A*δ sparse(δ*I, n, n) spzeros(n, n);
                      spzeros(n, 2*n) sparse(δ*I, n, n);
                      spzeros(n, 3*n)])
        else
            P = expmat(Matrix([A*δ sparse(δ*I, n, n) spzeros(n, n);
                           spzeros(n, 2*n) sparse(δ*I, n, n);
                           spzeros(n, 3*n)]))
        end
        Phi1Adelta = P[1:n, (n+1):2*n]
    end

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
    discr_bloat_interpolation(cont_sys, δ, approx_model, pade_expm, lazy_expm)

Compute bloating factors using forward or backward interpolation.

## Input

- `cs`           -- a continuous system
- `δ`            -- step size
- `approx_model` -- choose the approximation model among `"forward"` and
                    `"backward"`
- `pade_expm`    -- if true, use Pade approximant method to compute the
                    matrix exponential
- `lazy_expm`   --  if true, compute the matrix exponential in a lazy way
                    suitable for very large systems)

## Algorithm

See Frehse et al., CAV'11, *SpaceEx: Scalable Verification of Hybrid Systems*,
Lemma 3.

Note that in the unlikely case that A is invertible, the result can also
be obtained directly, as a function of the inverse of A and `e^{At} - I`.

The matrix `P` is such that: `ϕAabs = P[1:n, 1:n]`,
`Phi1Aabsdelta = P[1:n, (n+1):2*n]`, and `Phi2Aabs = P[1:n, (2*n+1):3*n]`.
"""
function discr_bloat_interpolation(cont_sys::InitialValueProblem{<:AbstractContinuousSystem},
                                   δ::Float64,
                                   approx_model::String="forward",
                                   pade_expm::Bool=false,
                                   lazy_expm::Bool=false,
                                   lazy_sih::Bool=true)

    sih = lazy_sih ? SymmetricIntervalHull : symmetric_interval_hull

    # unrwap coefficients matrix and initial states
    A, X0 = cont_sys.s.A, cont_sys.x0
    n = size(A, 1)

    # compute matrix ϕ = exp(Aδ)
    ϕ = exp_Aδ(A, δ, lazy_expm, pade_expm)

    # early return for homogeneous systems
    if cont_sys isa IVP{<:LinearContinuousSystem}
         Ω0 = CH(X0, ϕ * X0)
        return DiscreteSystem(ϕ, Ω0)
    end
    U = inputset(cont_sys)
    inputs = next_set(U, 1)

    # compute the transformation matrix to bloat the initial states
    if lazy_expm
        P = SparseMatrixExp([abs.(A*δ) sparse(δ*I, n, n) spzeros(n, n);
                             spzeros(n, 2*n) sparse(δ*I, n, n);
                             spzeros(n, 3*n)])
        Phi2Aabs = sparse(get_columns(P, (2*n+1):3*n)[1:n, :])
    else
        if pade_expm
            P = padm([abs.(A*δ) sparse(δ*I, n, n) spzeros(n, n);
                      spzeros(n, 2*n) sparse(δ*I, n, n);
                      spzeros(n, 3*n)])
        else
            P = expmat(Matrix([abs.(A*δ) sparse(δ*I, n, n) spzeros(n, n);
                           spzeros(n, 2*n) sparse(δ*I, n, n);
                           spzeros(n, 3*n)]))
        end
        Phi2Aabs = P[1:n, (2*n+1):3*n]
    end

    if isa(inputs, ZeroSet)
        if approx_model == "forward" || approx_model == "backward"
            Ω0 = CH(X0, ϕ * X0 + δ * inputs)
        end
    else
        EPsi = sih(Phi2Aabs * sih(A * inputs))
        discretized_U = δ * inputs + EPsi
        if approx_model == "forward"
            EOmegaPlus = sih(Phi2Aabs * sih((A * A) * X0))
            Ω0 = CH(X0, ϕ * X0 + discretized_U + EOmegaPlus)
        elseif approx_model == "backward"
            EOmegaMinus = sih(Phi2Aabs * sih((A * A * ϕ) * X0))
            Ω0 = CH(X0, ϕ * X0 + discretized_U + EOmegaMinus)
        end
    end

    if U isa ConstantInput
        return DiscreteSystem(ϕ, Ω0, discretized_U)
    else
        discretized_U = [δ * Ui + sih(Phi2Aabs * sih(A * Ui)) for Ui in U]
        return DiscreteSystem(ϕ, Ω0, discretized_U)
    end
end

function discr_bloat_interpolation(cont_sys::InitialValueProblem{<:LinearContinuousSystem},
                                   δ::Float64,
                                   approx_model::String="forward",
                                   pade_expm::Bool=false,
                                   lazy_expm::Bool=false,
                                   lazy_sih::Bool=true)

    # unwrap coefficients matrix and initial states
    A, X0 = cont_sys.s.A, cont_sys.x0

    # compute matrix ϕ = exp(Aδ)
    ϕ = exp_Aδ(A, δ, lazy_expm, pade_expm)

    
end
