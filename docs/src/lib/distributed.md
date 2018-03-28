# Distributed computations

This section of the manual describes functions to make use of distributed computation.

```@contents
Pages = ["distributed.md"]
```

```@meta
CurrentModule = Reachability
```

## Using multiple threads

To control the number of threads used by your BLAS library, use the function
`Base.LinAlg.BLAS.set_num_threads(n)`, where `n` is an integer. Furthermore,
the function `get_num_threads()` defined below will return the current value.

```julia
#
# This function is a part of Julia. License is MIT: https://julialang.org/license
#
function get_num_threads() # anonymous so it will be serialized when called
    blas = LinearAlgebra.BLAS.vendor()
    # Wrap in a try to catch unsupported blas versions
    try
        if blas == :openblas
            return ccall((:openblas_get_num_threads, Base.libblas_name), Cint, ())
        elseif blas == :openblas64
            return ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
        elseif blas == :mkl
            return ccall((:MKL_Get_Max_Num_Threads, Base.libblas_name), Cint, ())
        end

        # OSX BLAS looks at an environment variable
        if Sys.isapple()
            return ENV["VECLIB_MAXIMUM_THREADS"]
        end
    end

    return nothing
end
```
