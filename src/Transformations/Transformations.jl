__precompile__()
"""
Module to apply coordinate transformations.
"""
module Transformations

using LazySets, ..Systems

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
function transform(S::T;
                   method::String="schur")::Tuple{T, AbstractMatrix} where
                                            {T<:Union{DiscreteSystem, ContinuousSystem}}

    if method == "schur"
        return schur_transform(S)
    else
        error("The transformation method $method is undefined")
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
function schur_transform(S::T)::Tuple{T, AbstractMatrix} where
                                     {T<:Union{DiscreteSystem, ContinuousSystem}}

    A_new, T_new = schur(full(S.A)) # full (dense) matrix is required


    # recall that for Schur matrices, inv(T) == T'
    Z_inverse = T_new.'

    # apply transformation to the initial states
    X0_new = Z_inverse * S.X0

    # apply transformation to the inputs
    U_new = Z_inverse * S.U

    # obtain the transformation matrix for reverting the transformation again
    T_inverse = F[:vectors]

    if S isa DiscreteSystem
        return (DiscreteSystem(A_new, X0_new, S.δ, U_new), T_inverse)
    elseif S isa ContinuousSystem
        return (ContinuousSystem(A_new, X0_new, U_new), T_inverse)
    end
end

end # module
