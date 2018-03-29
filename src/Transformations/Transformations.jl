__precompile__()
"""
Module to apply coordinate transformations.
"""
module Transformations

using LazySets, ..Utils, Systems

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

    A_new, T_new = schur(full(S.A)) # full (dense) matrix is required


    # recall that for Schur matrices, inv(T) == T'
    Z_inverse = T_new.'

    # apply transformation to the initial states
    X0_new = Z_inverse * S.x0

    # apply transformation to the inputs
    U_new = Z_inverse * S.s.U

    # obtain the transformation matrix for reverting the transformation again
    T_inverse = F[:vectors]

    if S.s isa DiscreteSystem
        return (DiscreteSystem(A_new, X0_new, U_new), T_inverse)
    elseif S.s isa ContinuousSystem
        return (ContinuousSystem(A_new, X0_new, U_new), T_inverse)
    end
end

end # module
