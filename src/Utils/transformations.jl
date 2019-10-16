using LinearAlgebra

export transform,
       backtransform

"""
    backtransform(Rsets::ReachSolution, options::Options)

Undo a coordinate transformation.

### Input

- `Rsets`  -- flowpipe
- `option` -- problem options containing an `:transformation_matrix` entry

### Output

A new flowpipe where each reach set has been transformed.

### Notes

The transformation is implemented with a lazy `LinearMap`.
"""
function backtransform(Rsets, options::Options)
    transformation_matrix = options[:transformation_matrix]
    if transformation_matrix == nothing
        return Rsets
    end
    return project(Rsets, transformation_matrix)
end

"""
    transform(problem::InitialValueProblem, options::Options)

Interface function that calls the respective transformation function.

### Input

- `problem` -- discrete or continuous initial-value problem
- `option`  -- problem options

### Output

A tuple containing the transformed problem and the transformed options.

### Notes

The functions that are called in the background should return a the transformed
system components `A`, `X0`, and `U`, and also an inverse transformation matrix `M`.
If the system has an invariant, it is transformed as well.
"""
function transform(problem::InitialValueProblem, options::Options)
    method = options[:coordinate_transformation]
    if method == ""
        nothing # no-op
    elseif method == "schur"
        problem, T_inverse = schur_transform(problem)
        options[:transformation_matrix] = T_inverse
    else
        error("the transformation method $method is undefined")
    end
    return (problem, options)
end

"""
    schur_transform(problem::InitialValueProblem)

Applies a Schur transformation to a discrete or continuous initial-value problem.

### Input

- `problem` -- discrete or continuous initial-value problem

### Output

Transformed problem.

### Algorithm

We use Julia's default `schurfact` function to compute a
[Schur decomposition](https://en.wikipedia.org/wiki/Schur_decomposition)
of the coefficients matrix ``A``.
"""
function schur_transform(problem::InitialValueProblem{PT, ST}
                         ) where {PT <: Union{ConstrainedLinearControlDiscreteSystem, ConstrainedLinearControlContinuousSystem}, ST<:LazySet}

    n = size(problem.s.A, 1)

    # full (dense) matrix is required by schur
    # result S is a struct such that
    # - S.Schur is in Schur form and
    # - A == S.vectors * S.Schur * S.vectors'
    S = schur(Matrix(problem.s.A))

    # recall that for Schur matrices, inv(T) == T'
    Z_inverse = copy(transpose(S.vectors))

    # obtain transformed system matrix
    A_new = S.Schur

    # apply transformation to the initial states
    X0_new = Z_inverse * problem.x0

    # apply transformation to the inputs
    B_new = Z_inverse * problem.s.B
    U_new = problem.s.U

    # matrix for reverting the transformation again
    T_inverse = S.vectors

    # apply transformation to the invariant
    if hasmethod(stateset, Tuple{typeof(problem.s)})
        invariant_new = T_inverse * problem.s.X
    else
        invariant_new = Universe(n)
    end

    system_new = _wrap_system(PT, A_new, B_new, invariant_new, U_new)
    problem_new = InitialValueProblem(system_new, X0_new)
    return problem_new, T_inverse
end

function _wrap_system(PT::Type{<:ConstrainedLinearControlDiscreteSystem},
                      A, B, invariant, U)
    return ConstrainedLinearControlDiscreteSystem(A, B, invariant, U)
end

function _wrap_system(PT::Type{<:ConstrainedLinearControlContinuousSystem},
                      A, B, invariant, U)
    return ConstrainedLinearControlContinuousSystem(A, B, invariant, U)
end
