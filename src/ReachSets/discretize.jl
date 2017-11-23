export discretize

"""
    discretize(cont_sys, δ; [approx_model], [pade_expm], [lazy_expm])

Discretize a continuous system of ODEs with nondeterministic inputs.

## Input

- `cont_sys`          -- continuous system
- `δ`                 -- step size
- `approx_model`      -- the method to compute the approximation model for the
                         discretization, among:

    - `forward`    -- use forward-time interpolation
    - `backward`   -- use backward-time interpolation
    - `firstorder` -- use first order approximation of the ODE
    - `nobloating` -- do not bloat the initial states
                      (use for discrete-time reachability)

- `pade_expm`         -- (optional, default = `false`) if true, use Pade approximant
                         method to compute matrix exponentials of sparse matrices;
                         otherwise use Julia's buil-in `expm`
- `lazy_expm`         -- (optional, default = `false`) if true, compute the matrix
                         exponential in a lazy way (suitable for very large systems)

## Output

A discrete system.

## Notes

This function applies an approximation model to transform a continuous affine system
into a discrete affine system. This transformation allows to do dense time reachability,
i.e. such that the trajectories of the given continuous system are included in the
computed flowpipe of the discretized system.
For discrete-time reachability, use `approx_model="nobloating"`.
"""
function discretize(cont_sys::ContinuousSystem, δ::Float64;
                    approx_model::String="forward",
                    pade_expm::Bool=false,
                    lazy_expm::Bool=false)::DiscreteSystem

    if approx_model in ["forward", "backward"]
        return discr_bloat_interpolation(cont_sys, δ, approx_model, pade_expm, lazy_expm)
    elseif approx_model == "firstorder"
        return discr_bloat_firstorder(cont_sys, δ)
    elseif approx_model == "nobloating"
        return discr_no_bloat(cont_sys, δ, pade_expm, lazy_expm)
    else
        error("The approximation model is invalid")
    end
end

"""
    bloat_firstorder(cont_sys, δ)

Compute bloating factors using first order approximation.

## Input

- `cont_sys` -- a continuous affine system
- `δ`        -- step size

## Notes

In this algorithm, the infinity norm is used.
See also: `discr_bloat_interpolation` for more accurate (less conservative)
bounds.

## Algorithm

This uses a first order approximation of the ODE, and matrix norm upper bounds,
see Le Guernic, C., & Girard, A., 2010, *Reachability analysis of linear systems
using support functions. Nonlinear Analysis: Hybrid Systems, 4(2), 250-262.*
"""
function discr_bloat_firstorder(cont_sys::ContinuousSystem, δ::Float64)::DiscreteSystem

    if !(cont_sys.U isa ConstantNonDeterministicInput)
        error("This discretization algorithm is only implemented for constant inputs")
    end

    Anorm = norm(full(cont_sys.A), Inf)
    RX0 = norm(cont_sys.X0, Inf)
    inputs = get_set(cont_sys.U)
    RU = norm(inputs, Inf)
    α = (exp(δ*Anorm) - 1. - δ*Anorm)*(RX0 + RU/Anorm)
    β = (exp(δ*Anorm) - 1. - δ*Anorm)*RU/Anorm
    ϕ = expm(full(cont_sys.A))
    Ω0 = CH(cont_sys.X0, ϕ * cont_sys.X0 + δ*inputs + Ball2(zeros(size(ϕ, 1)), α))
    discr_U =  δ * inputs + Ball2(zeros(size(ϕ, 1)), β)
    return DiscreteSystem(ϕ, Ω0, δ, discr_U)
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
function discr_no_bloat(cont_sys::ContinuousSystem,
                        δ::Float64,
                        pade_expm::Bool,
                        lazy_expm::Bool)::DiscreteSystem

    n = size(cont_sys.A, 1)
    if lazy_expm
        ϕ = SparseMatrixExp(cont_sys.A * δ)
    else
        if pade_expm
            ϕ = padm(cont_sys.A * δ)
        else
            ϕ = expm(full(cont_sys.A * δ))
        end
    end

    # early return for homogeneous systems
    inputs = get_set(cont_sys.U, 1)
    if isa(inputs, VoidSet) && length(cont_sys.U) == 1
            Ω0 = cont_sys.X0
            return DiscreteSystem(ϕ, Ω0, δ)
    end

    # compute matrix to transform the inputs
    if lazy_expm
        P = SparseMatrixExp([cont_sys.A*δ sparse(δ*I, n, n) spzeros(n, n); spzeros(n, 2*n) sparse(δ*I, n, n); spzeros(n, 3*n)])
        Phi1Adelta = get_columns(P, (n+1):2*n)[1:n, :]
    else
        if pade_expm
            P = padm([cont_sys.A*δ sparse(δ*I, n, n) spzeros(n, n); spzeros(n, 2*n) sparse(δ*I, n, n); spzeros(n, 3*n)])
        else
            P = expm(full([cont_sys.A*δ sparse(δ*I, n, n) spzeros(n, n); spzeros(n, 2*n) sparse(δ*I, n, n); spzeros(n, 3*n)]))
        end
        Phi1Adelta = P[1:n, (n+1):2*n]
    end

    discretized_U = Phi1Adelta * inputs

    Ω0 = cont_sys.X0

    if length(cont_sys.U) == 1
        return DiscreteSystem(ϕ, Ω0, δ, discretized_U)
    else
        discretized_U_arr = Vector{LazySet}(length(cont_sys.U))
        discretized_U_arr[1] = discretized_U
        for i in 2:length(cont_sys.U)
            inputs = next(cont_sys.U, i)[1]
            discretized_U_arr[i] = Phi1Adelta * inputs
        end
        return DiscreteSystem(ϕ, Ω0, δ, discretized_U_arr)
    end
end

"""
    discr_bloat_interpolation(cont_sys, δ, approx_model, pade_expm, lazy_expm)

Compute bloating factors using forward or backward interpolation.

## Input

- `cs`           -- a continuous system
- `δ`            -- step size
- `approx_model` -- choose the approximation model among `"forward"` and `"backward"`
- `pade_expm`    -- if true, use Pade approximant method to compute the
                    matrix exponential
- `lazy_expm`   --  if true, compute the matrix exponential in a lazy way
                    suitable for very large systems)

## Algorithm

See Frehse et al CAV'11 paper, *SpaceEx: Scalable Verification of Hybrid Systems*,
see Lemma 3.

Note that in the unlikely case that A is invertible, the result can also
be obtained directly, as a function of the inverse of A and `e^{At} - I`.

The matrix `P` is such that: `ϕAabs = P[1:n, 1:n]`, `Phi1Aabsdelta = P[1:n, (n+1):2*n]`,
and `Phi2Aabs = P[1:n, (2*n+1):3*n]`.
"""
function discr_bloat_interpolation(cont_sys::ContinuousSystem,
                                   δ::Float64,
                                   approx_model::String,
                                   pade_expm::Bool,
                                   lazy_expm::Bool)::DiscreteSystem

    n = size(cont_sys.A, 1)

    # compute matrix ϕ = exp(Aδ)
    if lazy_expm
        ϕ = SparseMatrixExp(cont_sys.A*δ)
    else
        if pade_expm
            ϕ = padm(cont_sys.A*δ)
        else
            ϕ = expm(full(cont_sys.A*δ))
        end
    end

    # early return for homogeneous systems
    inputs = get_set(cont_sys.U, 1)
    if isa(inputs, VoidSet) && length(cont_sys.U) == 1
            Ω0 = CH(cont_sys.X0, ϕ * cont_sys.X0)
            return DiscreteSystem(ϕ, Ω0, δ)
    end

    # compute the transformation matrix to bloat the initial states
    if lazy_expm
        P = SparseMatrixExp([abs.(cont_sys.A*δ) sparse(δ*I, n, n) spzeros(n, n); spzeros(n, 2*n) sparse(δ*I, n, n); spzeros(n, 3*n)])
        Phi2Aabs = get_columns(P, (2*n+1):3*n)[1:n, :]
    else
        if pade_expm
            P = padm([abs.(cont_sys.A*δ) sparse(δ*I, n, n) spzeros(n, n); spzeros(n, 2*n) sparse(δ*I, n, n); spzeros(n, 3*n)])
        else
            P = expm(full([abs.(cont_sys.A*δ) sparse(δ*I, n, n) spzeros(n, n); spzeros(n, 2*n) sparse(δ*I, n, n); spzeros(n, 3*n)]))
        end
        Phi2Aabs = P[1:n, (2*n+1):3*n]
    end

    if isa(inputs, VoidSet)
        if approx_model == "forward"
            Ω0 = CH(cont_sys.X0, ϕ * cont_sys.X0 + δ * inputs)
        elseif approx_model == "backward"
            Ω0 = CH(cont_sys.X0, ϕ * cont_sys.X0 + δ * inputs)
        end
    else
        EPsi = symmetric_interval_hull(Phi2Aabs * symmetric_interval_hull(cont_sys.A * inputs))
        discretized_U = δ * inputs + EPsi
        if approx_model == "forward"
            EOmegaPlus = symmetric_interval_hull(Phi2Aabs * symmetric_interval_hull((cont_sys.A * cont_sys.A) * cont_sys.X0))
            Ω0 = CH(cont_sys.X0, ϕ * cont_sys.X0 + discretized_U + EOmegaPlus)
        elseif approx_model == "backward"
            EOmegaMinus = symmetric_interval_hull(Phi2Aabs * symmetric_interval_hull((cont_sys.A * cont_sys.A * ϕ) * cont_sys.X0))
            Ω0 = CH(cont_sys.X0, ϕ * cont_sys.X0 + discretized_U + EOmegaMinus)
        end
    end

    if length(cont_sys.U) == 1
        return DiscreteSystem(ϕ, Ω0, δ, discretized_U)
    else
        discretized_U_arr = Vector{LazySet}(length(cont_sys.U))
        discretized_U_arr[1] = discretized_U
        for i in 2:length(cont_sys.U)
            inputs = next(cont_sys.U, i)[1]
            EPsi_i = symmetric_interval_hull(Phi2Aabs * symmetric_interval_hull(cont_sys.A * inputs))
            discretized_U_arr[i] = δ * inputs + EPsi_i
        end
        return DiscreteSystem(ϕ, Ω0, δ, discretized_U_arr)
    end
end
