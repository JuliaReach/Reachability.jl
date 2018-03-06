"""
    check_property(S, N; [algorithm], [ε_init], [set_type_init], [ε_iter],
                   [set_type_iter], [assume_sparse], [assume_homogeneous],
                   [lazy_X0], [kwargs]...)

Interface to property checking algorithms for an LTI system.

### Input

- `S`                  -- LTI system, discrete or continuous
- `N`                  -- number of computed sets
- `algorithm`          -- (optional, default: `"explicit"`), algorithm backend;
                          see `available_algorithms` for all admissible options
- `ε_init`             -- (optional, default: `Inf`) error bound for the
                          approximation of the initial states (during
                          decomposition)
- `set_type_init`      -- (optional, default: `Hyperrectangle`) set type for the
                          approximation of the initial states (during
                          decomposition)
- `ε_iter`             -- (optional, default: `Inf`) error bound for the
                          approximation of the inputs
- `set_type_iter`      -- (optional, default: `Hyperrectangle`) set type for the
                          approximation of the inputs
- `assume_sparse`      -- (optional, default: `true`) if true, it is assumed
                          that the coefficients matrix (exponential) is sparse;
                          otherwise, it is transformed to a full matrix
- `assume_homogeneous` -- (optional, default: `false`) if true, it is assumed
                          that the system has no input (linear system), and it
                          case it has one, the input is ignored
- `lazy_X0`            -- (optional, default: `false`) if true, transform the
                          set of initial states to the caretsian product of
                          two-dimensional polygons; otherwise, the given input,
                          as a lazy set, is passed to the backend
- `kwargs`             -- (optional) additional arguments that are passed to
                          the backend

### Notes

A dictionary with available algorithms is available via
`Properties.available_algorithms`.

WARNING:

Only systems of even dimension are parsed; for odd dimension, manually add an
extra variable with no dynamics.
"""
function check_property(S::AbstractSystem,
                        N::Int;
                        algorithm::String="explicit",
                        ε_init::Float64=Inf,
                        set_type_init::Type{<:LazySet}=Hyperrectangle,
                        ε_iter::Float64=Inf,
                        set_type_iter::Type{<:LazySet}=Hyperrectangle,
                        assume_sparse=true,
                        assume_homogeneous=false,
                        lazy_X0=false,
                        kwargs...)::Int

    # unpack arguments
    kwargs_dict = Dict(kwargs)

    # list containing the arguments passed to any reachability function
    args = []

    #coefficients matrix
    push!(args, S.A)

    # Cartesian decomposition of the initial set
    if lazy_X0
        Xhat0 = S.X0
    elseif !isempty(kwargs_dict[:block_types])
        Xhat0 = array(decompose(S.X0, ɛ=ε_init,
                                block_types=kwargs_dict[:block_types]))
    elseif set_type_init == LazySets.Interval
        Xhat0 = array(decompose(S.X0, set_type=set_type_init, ɛ=ε_init,
                                blocks=ones(Int, dim(S.X0))))
    else
        Xhat0 = array(decompose(S.X0, set_type=set_type_init, ɛ=ε_init))
    end

    # shortcut if only the initial set is required
    if N == 1
        if length(kwargs_dict[:blocks]) == 1
            Xhat0_mod = Xhat0[kwargs_dict[:blocks][1]]
        else
            Xhat0_mod = CartesianProductArray(Xhat0)
        end
        return check_property(Xhat0_mod, kwargs_dict[:property]) ? 0 : 1
    end
    push!(args, Xhat0)

    if !assume_homogeneous
        push!(args, S.U)

        # overapproximation function (with or without iterative refinement)
        if ε_iter < Inf
            push!(args, x -> overapproximate(x, set_type_iter, ε_iter))
        else
            push!(args, x -> overapproximate(x, set_type_iter))
        end
    end

    # ambient dimension
    n = Systems.dim(S)
    push!(args, n)

    # size of each block
    @assert (n % 2 == 0) "the number of dimensions should be even"
    push!(args, div(n, 2))

    # number of computed sets
    push!(args, N)

    # add mode-specific block(s) argument
    if algorithm == "explicit"
        if !(:blocks in keys(kwargs_dict))
            error("This algorithm needs a specified block argument.")
        elseif length(kwargs_dict[:blocks]) == 1
            bi = kwargs_dict[:blocks][1]
            push!(args, bi)
            algorithm_backend = "explicit_block"
        else
            push!(args, kwargs_dict[:blocks])
            push!(args, kwargs_dict[:partition])
            algorithm_backend = "explicit_blocks"
        end
    else
        error("Unsupported algorithm: ", algorithm)
    end

    # add property
    push!(args, kwargs_dict[:property])

    # call the adequate function with the given arguments list
    answer = available_algorithms[algorithm_backend]["func"](args...)

    # return the result
    return answer
end
