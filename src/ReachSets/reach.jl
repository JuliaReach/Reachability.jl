"""
    reach(S, N; [algorithm], [ɛ], [assume_sparse], [assume_homogeneous],
          [set_type], [numeric_type], [lazy_X0], [kwargs]...)

Interface to reachability algorithms for an LTI system.

### Input

- `S`                  -- LTI system, discrete or continuous
- `N`                  -- number of computed sets
- `algorithm`          -- (optional, default: `"explicit"`), reachability
                          algorithm backend; see `available_algorithms` for all
                          admissible options
- `ɛ`                  -- (optional, default: `Inf`) vector for error tolerance
                          on each block
- `assume_sparse`      -- (optional, default: `true`) if true, it is assumed
                          that the coefficients matrix (exponential) is sparse;
                          otherwise, it is transformed to a full matrix
- `assume_homogeneous` -- (optional, default: `false`) if true, it is assumed
                          that the system has no input (linear system), and in
                          case it has one, the input is ignored
- `set_type`           -- (optional, default: `Hyperrectangle`) type of set that
                          is used for overapproximation in block dimensions
- `numeric_type`       -- (optional, default: `Float64`) numeric type of the
                          resulting set
- `lazy_X0`            -- (optional, default: `false`) if true, transform the
                          set of initial states to the caretsian product of
                          two-dimensional polygons; otherwise, the given input,
                          as a lazy set, is passed to the backend
- `kwargs`             -- (optional) additional arguments that are passed to the
                          backend

### Notes

A dictionary with available algorithms is available via `Reachability.available_algorithms`.

WARNING:

Only systems of even dimension are parsed; for odd dimension, manually add an
extra variable with no dynamics.
"""
function reach(S::AbstractSystem,
               N::Int;
               algorithm::String="explicit",
               ɛ::Float64=Inf,
               assume_sparse=true,
               assume_homogeneous=false,
               set_type::Type=Hyperrectangle,
               numeric_type::Type=Float64,
               lazy_X0=false,
               kwargs...)::Union{Vector{<:CartesianProductArray}, Vector{<:set_type}}

    # unpack arguments
    kwargs_dict = Dict(kwargs)

    # list containing the arguments passed to any reachability function
    args = []

    #coefficients matrix
    if assume_sparse
        push!(args, sparse(S.A))
    else
        push!(args, S.A)
    end

    # Cartesian decomposition of the initial set
    if lazy_X0
        Xhat0 = S.X0
    else
        Xhat0 = array(decompose(S.X0, set_type=set_type, ɛ=ɛ))
    end

    # shortcut if only the initial set is required
    if N == 1
        if length(kwargs_dict[:blocks]) == 1
            return [Xhat0[kwargs_dict[:blocks][1]]]
        else
            return Xhat0
        end
    end
    push!(args, Xhat0)

    if !assume_homogeneous
        push!(args, S.U)
    end

    # ambient dimension
    n = Systems.dim(S)
    push!(args, n)

    # size of each block
    @assert (n % 2 == 0) "the number of dimensions should be even"
    push!(args, div(n, 2))

    # number of computed sets
    push!(args, N)

    # overapproximation function (with or without iterative refinement)
    if ɛ == Inf
        push!(args, x -> overapproximate(x, set_type))
    else
        push!(args, x -> overapproximate(x, ɛ))
    end

    # preallocate output vector and add mode-specific block(s) argument
    if algorithm == "explicit"
        if !(:blocks in keys(kwargs_dict))
            error("This algorithm needs a specified block argument.")
        elseif length(kwargs_dict[:blocks]) == 1
            bi = kwargs_dict[:blocks][1]
            push!(args, bi)
            res = Vector{set_type{numeric_type}}(N)
            algorithm_backend = "explicit_block"
        else
            blocks = kwargs_dict[:blocks]
            push!(args, blocks)
            res = Vector{CartesianProductArray{numeric_type}}(N)
            algorithm_backend = "explicit_blocks"
        end
    else
        error("Unsupported algorithm: ", algorithm)
    end
    push!(args, res)

    # call the adequate function with the given arguments list
    available_algorithms[algorithm_backend]["func"](args...)

    # return the result
    return res
end
