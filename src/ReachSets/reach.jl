"""
    reach(S, N; [algorithm], [ɛ], [iterative_refinement], [assume_sparse],
                [assume_homogeneous], [kwargs]...)

Interface to reachability algorithms for an affine system with non-deterministic inputs.

INPUT:

- `S`                      -- affine system, discrete or continuous
- `N`                      -- number of computed sets
- `algorithm`              -- (optional, default: `explicit`), reachability
                              algorithm backend; see `available_algorithms` for all
                              admissible options
- `ɛ`                      -- (optional, default: Inf) vector for error tolerances on each
                              block; if iterative_refinement is set to false, this
                              value is only used to decompose the initial states,
                              and ignored afterwards
- `iterative_refinement`   -- (optional default: `false`) if true, perform iterative
                              refinement with the given tolerance ɛ; otherwise,
                              only box directions are computed
- `assume_sparse`          -- (optional, default: `true`) if true, it is assumed that the
                              coefficients matrix (exponential) is sparse; otherwise,
                              it is transformed to a full matrix
- `assume_homogeneous`     -- (optional, default: `false`) if true, it is assumed that the
                              system has no input (linear system), and it case it has one,
                              the input is ignored; otherwise, the given non-deterministic
                              input is passed to the backend
- `lazy_X0`                -- (optional, default: `false`) if true, transform the
                              set of initial states to the caretsian product of
                              two-dimensional polygons; otherwise, the given input,
                              as a lazy set, is passed to the backend 
- `kwargs`                 -- (optional) additional arguments that are passed to the backend

NOTES:

- A dictionary with available algorithms is available at `Reachability.available_algorithms`.

WARNING:

- Only systems of even dimension are parsed; if it is not your case, manually
  add an extra variable with no dynamics.
"""
function reach(S::Union{DiscreteSystem, ContinuousSystem},
               N::Int64;
               algorithm::String="explicit",
               ɛ::Float64=Inf,
               iterative_refinement=false,
               assume_sparse=true,
               assume_homogeneous=false,
               lazy_X0=false,
               kwargs...)::Union{Vector{CartesianProductArray}, Vector{HPolygon}}

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
       if iterative_refinement
            Xhat0 = decompose(S.X0, ɛ)
        else
            Xhat0 = decompose(S.X0)
        end
        Xhat0 = convert(Vector{LazySets.HPolygon}, Xhat0.sfarray)
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

    # number of blocks
    b = div(n, 2)
    if n % 2 != 0
        error("the number of dimensions should be even")
    end
    push!(args, b)

    # number of computed sets
    push!(args, N)

    # overapproximation function (with or without iterative refinement)
    if iterative_refinement
        push!(args, x -> overapproximate(x, ɛ))
    else
        push!(args, x -> overapproximate(x))
    end

    # preallocate output vector and add mode-specific block(s) argument
    if algorithm == "explicit"
        if !(:blocks in keys(kwargs_dict))
            error("This algorithm needs a specified block argument.")
        elseif length(kwargs_dict[:blocks]) == 1
            bi = kwargs_dict[:blocks][1]
            push!(args, bi)
            res = Vector{HPolygon}(N)
            algorithm_backend = "explicit_block"
        else
            blocks = kwargs_dict[:blocks]
            push!(args, blocks)
            res = Vector{CartesianProductArray}(N)
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