"""
    check_property(S, N; [algorithm], [approx_init], [approx_sets],
                   [assume_sparse], [assume_homogeneous], [lazy_X0],
                   [kwargs]...)

Interface to property checking algorithms for an LTI system.

### Input

- `S`                  -- LTI system, discrete or continuous
- `N`                  -- number of computed sets
- `algorithm`          -- (optional, default: `"explicit"`), algorithm backend;
                          see `available_algorithms` for all admissible options
- `approx_init`        -- (optional, default: `Hyperrectangle`) set
                          approximation parameter for the initial states (during
                          decomposition)
- `approx_sets`        -- (optional, default: `Hyperrectangle`) set
                          approximation parameter for the inputs
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

A dictionary with available algorithms is available via `Properties.available_algorithms`.

WARNING:

Only systems of even dimension are parsed; for odd dimension, manually add an
extra variable with no dynamics.
"""
function check_property(S::AbstractSystem,
                        N::Int;
                        algorithm::String="explicit",
                        approx_init::Union{Float64, Type{<:LazySet}}=Hyperrectangle,
                        approx_sets::Union{Float64, Type{<:LazySet}}=Hyperrectangle,
                        assume_sparse=true,
                        assume_homogeneous=false,
                        lazy_X0=false,
                        kwargs...)::Int

    # unpack arguments
    kwargs_dict = Dict(kwargs)

    # list containing the arguments passed to any reachability function
    args = []

    #coefficients matrix
    if assume_sparse
        push!(args, S.A)
    else
        try
            push!(args, full(S.A))
        catch
            push!(args, S.A)
        end
    end

    # Cartesian decomposition of the initial set
    if lazy_X0
        Xhat0 = S.X0
    elseif approx_init isa Float64
        Xhat0 = array(decompose(S.X0, ɛ=approx_init))
    else
        Xhat0 = array(decompose(S.X0, set_type=approx_init))
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
        if approx_sets isa Float64
            push!(args, x -> overapproximate(x, HPolygon, approx_sets))
        else
            push!(args, x -> overapproximate(x, approx_sets))
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
            blocks = kwargs_dict[:blocks]
            push!(args, blocks)
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
