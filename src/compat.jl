#=
This file defines internal functions for compatibility across
different Julia versions.
=#

using Compat
using Compat: copyto!, axes, argmax
import Compat.String
using Compat.LinearAlgebra
import Compat.LinearAlgebra: norm, checksquare, LAPACKException,
                             SingularException, eye, ×
import Compat.InteractiveUtils.subtypes
export _At_mul_B, _A_mul_B!

@static if VERSION < v"0.7-"
    @inline _At_mul_B(A, B) = At_mul_B(A, B)
    @inline _A_mul_B!(C, A, B) = A_mul_B!(C, A, B)
    expmat = expm
else
    using SparseArrays, Printf

    @inline _At_mul_B(A, B) = transpose(A) * B
    @inline _A_mul_B!(C, A, B) = mul!(C, A, B)
    expmat = exp
end

if VERSION > v"1.0-"
    export eye
end
