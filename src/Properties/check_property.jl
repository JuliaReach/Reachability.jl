"""
    check_property(S, N; [algorithm], [ɛ], [iterative_refinement], [assume_sparse],
                [assume_homogeneous], [kwargs]...)

Property checking algorithm of an affine system.

INPUT:

- ``S``                    -- affine system, discrete or continuous
- ``N``                    -- number of computed sets
- ``algorithm``            -- (optional, default: ``'explicit'``), algorithm backend;
                              see ``available_algorithms`` for all admissible options
- ``ɛ``                    -- (optional, default: Inf) vector for error tolerances on
                              each block; if iterative_refinement is set to false,
                              this value is only used to decompose the initial states,
                              and ignored afterwards
- ``iterative_refinement`` -- (optional default: false) if true, perform iterative
                              refinement with the given tolerance ɛ; otherwise,
                              only box directions are computed; currently we do not
                              support this option
- ``assume_sparse``        -- (optional, default: true) if true, it is assumed that the
                              coefficients matrix (exponential) is sparse; otherwise,
                              it is transformed to a full matrix
- ``assume_homogeneous``   -- (optional, default: false) if true, it is assumed that the
                              system has no input (linear system), and it case it has one,
                              the input is ignored; otherwise, the given non-deterministic
                              input is passed to the backend
- `set_type`               -- (optional, default: `HPolygon`) type of set that is used
                              for overapproximation in 2D
- `lazy_X0`                -- (optional, default: `false`) if true, transform the
                              set of initial states to the caretsian product of
                              two-dimensional polygons; otherwise, the given input,
                              as a lazy set, is passed to the backend
- ``kwargs``               -- (optional) additional arguments that are passed to the backend

NOTES:

- A dictionary with available algorithms is available at ``Properties.available_algorithms``.

WARNING:

- Only systems of even dimension are parsed; for odd dimension please manually
  add an extra variable with no dynamics.
"""
function check_property(S::Union{DiscreteSystem, ContinuousSystem},
                        N::Int;
                        algorithm::String="explicit",
                        ɛ::Float64=Inf,
                        iterative_refinement=false,
                        assume_sparse=true,
                        assume_homogeneous=false,
                        set_type::Type=HPolygon,
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
    else
        if iterative_refinement
            Xhat0 = decompose(S.X0, ɛ, set_type)
        else
            Xhat0 = decompose(S.X0, set_type)
        end
        Xhat0 = array(Xhat0)
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
        if iterative_refinement
            push!(args, x -> overapproximate(x, ɛ))
        else
            push!(args, x -> overapproximate(x, set_type))
        end
    end

    # ambient dimension
    n = Systems.dim(S)
    push!(args, n)

    # number of blocks, b
    b = div(n, 2)
    if n % 2 != 0
        error("the number of dimensions should be even")
    end
    push!(args, b)

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
