__precompile__()

module Transformations

using LazySets, ..Systems

"""
    transform(method, Δ, plot_vars)

Forward call to the respective transformation function.

INPUT:

- ``method``    -- transformation method name, one of: ``'schur'``
- ``Δ``         -- discretized system
- ``plot_vars`` -- variables to plot; two-dimensional index array
"""
function transform(method::String, Δ::DiscreteSystem, plot_vars::Array{Int64, 1})::Tuple{DiscreteSystem, SparseMatrixCSC{Float64,Int64}}
    if method == "schur"
        return transformSchur(Δ, plot_vars)
    else
        error("undefined transformation")
    end
end

"""
    transformSchur(Δ, plot_vars)

Applies a Schur transformation to a discretized system Δ.
The output is again a discretized system together with a transformation matrix
for undoing the transformation.
The argument plot_vars is used to determine whether we will add another
dimension for time later and hence whether the transformation matrix should have
another dimension for time.

INPUT:

- ``Δ``         -- discretized system
- ``plot_vars`` -- variables to plot; two-dimensional index array
"""
function transformSchur(Δ::DiscreteSystem, plot_vars::Array{Int64, 1})::Tuple{DiscreteSystem, SparseMatrixCSC{Float64,Int64}}
    A::SparseMatrixCSC{Float64, Int64} = Δ.A
    F = schurfact(full(A))
    A_new::SparseMatrixCSC{Float64, Int64} = sparse(F[:Schur])
    Z_inverse = sparse(F[:vectors].') # for Schur matrix: inv(F) == F'

    # apply transformation to initial states
    X0_new = Z_inverse * Δ.X0

    # apply transformation to inputs
    U_new = Z_inverse * Δ.U
    
    # compute the transformation matrix for undoing the transformation again
    transformationMatrix = sparse(F[:vectors])
#=
    got_time = (plot_vars[1] == 0)
    if got_time
        # add another dimension for time: block matrix [M 0; 0 1]
        transformationMatrix = sparse(cat([1, 2], transformationMatrix, [1]))
    end
=#
    return (DiscreteSystem(A_new, X0_new, Δ.δ, U_new), transformationMatrix)
end

export transform

end