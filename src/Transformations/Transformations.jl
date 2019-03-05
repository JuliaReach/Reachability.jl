__precompile__()
"""
Module to apply coordinate transformations.
"""
module Transformations

using LazySets, ..Utils, MathematicalSystems

include("../compat.jl")

export transform

"""
    transform(S; [method])

Interface function that calls the respective transformation function.

### Input

- `S`      -- discrete or continuous system
- `method` -- (optional, default: `'schur'`) transformation method name; valid
              otions are:

    * `'schur'`

### Output

A tuple containing:

- transformed discrete or continuous system
- inverse transformation matrix for reverting the transformation

### Notes

The functions that are called in the background should return a the transformed
system components `A`, `X0`, and `U`, and also an inverse transformation matrix `M`.
"""
function transform(S::InitialValueProblem; method::String="schur")

    if method == "schur"
        return schur_transform(S)
    else
        error("the transformation method $method is undefined")
    end
end

"""
    schur_transform(S)

Applies a Schur transformation to a discrete or continuous system.

### Input

- `S` -- discrete or continuous system

### Output

A tuple containing:

- transformed discrete or continuous system
- inverse transformation matrix for reverting the transformation

### Algorithm

We use Julia's default `schurfact` function to compute a
[Schur decomposition](https://en.wikipedia.org/wiki/Schur_decomposition)
of the coefficients matrix ``A``.
"""
function schur_transform(S::InitialValueProblem)

    A_new, T_new = schur(Matrix(S.A)) # full (dense) matrix is required


    # recall that for Schur matrices, inv(T) == T'
    Z_inverse = copy(transpose(T_new))

    # apply transformation to the initial states
    X0_new = Z_inverse * S.x0

    # apply transformation to the inputs
    U_new = Z_inverse * S.s.U

    # obtain the transformation matrix for reverting the transformation again
    T_inverse = F[:vectors]

    if S.s isa DiscreteSystem
        s = ConstrainedLinearControlDiscreteSystem(A_new, Matrix(1.0I, size(A_new)), nothing, U_new)
        p = InitialValueProblem(s, X0_new)
    elseif S.s isa ContinuousSystem
        s = ConstrainedLinearControlContinuousSystem(A_new, Matrix(1.0I, size(A_new)), nothing, U_new)
        p = InitialValueProblem(s, X0_new)
    end
    return (p, T_inverse)
end

end # module
