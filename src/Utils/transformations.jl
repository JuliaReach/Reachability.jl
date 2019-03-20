export transform

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
        options[:transformation_matrix] = Tm
        problem = schur_transform(problem)
    else
        error("the transformation method $method is undefined")
    end
    return (problem, options)
end

"""
    schur_transform(problem::InitialValueProblem)

Applies a Schur transformation to a discrete or continuous initial-value problem.

### Input

- `S` -- discrete or continuous initial-value problem

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
    A_new, T_new = schur(Matrix(problem.s.A))

    # recall that for Schur matrices, inv(T) == T'
    Z_inverse = copy(transpose(T_new))

    # apply transformation to the initial states
    X0_new = Z_inverse * problem.x0

    # apply transformation to the inputs
    U_new = Z_inverse * problem.s.U

    # obtain the transformation matrix for reverting the transformation again
    T_inverse = F[:vectors]

    # apply transformation to the invariant
    if hasmethod(stateset, Tuple{typeof(problem.s)})
        invariant_new = T_inverse * problem.s.X
    else
        invariant_new = Universe(n)
    end

    system_new = (A_new, I(n), invariant_new, U_new)
    return InitialValueProblem(PT(system_new), X0_new)
end
