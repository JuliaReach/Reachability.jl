__precompile__()
"""
Module to apply coordinate transformations.
"""
module Transformations

using LazySets, ..Systems

export transform

"""
    transform(S, method)

Interface function that calls the respective transformation function.

### Input

- `S`      -- discrete or continuous system
- `method` -- transformation method name, one of: `'schur'`

### Output

- transformed discrete or continuous system
- inverse transformation matrix for reverting the transformation

### Notes

The functions that are called in the background should return a four-tuple
`(A, X0, U, M)` consisting of the transformed system components `A`, `X0`, and
`U`, and an inverse transformation matrix `M`.
"""
function transform(S::T,
                   method::String)::Tuple{T, SparseMatrixCSC{Float64, Int64}} where
                       {T<:Union{DiscreteSystem, ContinuousSystem}}
    if method == "schur"
        (A, X0, U, M) = transform_schur(S)
    else
        error("undefined transformation")
    end

    if S isa DiscreteSystem
        return (DiscreteSystem(A, X0, S.δ, U), M)
    else
        return (ContinuousSystem(A, X0, U), M)
    end
end

"""
    transform_schur(S)

Applies a Schur transformation to a discrete or continuous system S.

### Input

- `S` -- discrete or continuous system

### Output

- transformed system matrix
- transformed initial states
- transformed nondeterministic inputs
- inverse transformation matrix for reverting the transformation

### Algorithm

We use the default `schurfact` function to compute a Schur decomposition.
"""
function transform_schur(S::Union{DiscreteSystem, ContinuousSystem})::Tuple{
                             SparseMatrixCSC{Float64, Int64},
                             LazySet,
                             NonDeterministicInput,
                             SparseMatrixCSC{Float64,Int64}}
    A::SparseMatrixCSC{Float64, Int64} = S.A
    F = schurfact(full(A))
    A_new::SparseMatrixCSC{Float64, Int64} = sparse(F[:Schur])
    Z_inverse = sparse(F[:vectors].') # for Schur matrix: inv(F) == F'

    # apply transformation to initial states
    X0_new = Z_inverse * S.X0

    # apply transformation to inputs
    U_new = Z_inverse * S.U
    
    # compute the transformation matrix for reverting the transformation again
    inverse = sparse(F[:vectors])

    return (A_new, X0_new, U_new, inverse)
end

end # module
