"""
    reach(S, N; [algorithm], [ε_init], [set_type_init], [ε_iter],
          [set_type_iter], [assume_sparse], [assume_homogeneous],
          [numeric_type], [lazy_X0], [kwargs]...)

Interface to reachability algorithms for an LTI system.

### Input

- `S`                  -- LTI system, discrete or continuous
- `N`                  -- number of computed sets
- `algorithm`          -- (optional, default: `"explicit"`), reachability
                          algorithm backend; see `available_algorithms` for all
                          admissible options
- `ε_init`             -- (optional, default: `Inf`) error bound for the
                          approximation of the initial states (during
                          decomposition)
- `set_type_init`      -- (optional, default: `Hyperrectangle`) set type for the
                          approximation of the initial states (during
                          decomposition)
- `ε_iter`             -- (optional, default: `Inf`) error bound for the
                          approximation of the states ``X_k``, ``k>0``
- `set_type_iter`      -- (optional, default: `Hyperrectangle`) set type for the
                          approximation of the states ``X_k``, ``k>0``
- `assume_sparse`      -- (optional, default: `true`) if true, it is assumed
                          that the coefficients matrix (exponential) is sparse;
                          otherwise, it is transformed to a full matrix
- `assume_homogeneous` -- (optional, default: `false`) if true, it is assumed
                          that the system has no input (linear system), and in
                          case it has one, the input is ignored
- `numeric_type`       -- (optional, default: `Float64`) numeric type of the
                          resulting set
- `lazy_X0`            -- (optional, default: `false`) if true, transform the
                          set of initial states to the caretsian product of
                          two-dimensional polygons; otherwise, the given input,
                          as a lazy set, is passed to the backend
- `kwargs`             -- (optional) additional arguments that are passed to the
                          backend

### Notes

A dictionary with available algorithms is available via
`Reachability.available_algorithms`.
"""
function reach(S::AbstractSystem,
               N::Int;
               algorithm::String="explicit",
               ε_init::Float64=Inf,
               set_type_init::Type{<:LazySet}=Hyperrectangle,
               ε_iter::Float64=Inf,
               set_type_iter::Type{<:LazySet}=Hyperrectangle,
               assume_sparse=true,
               assume_homogeneous=false,
               numeric_type::Type=Float64,
               lazy_X0=false,
               kwargs...)::Vector{<:LazySet}

    # unpack arguments
    kwargs_dict = Dict(kwargs)

    # list containing the arguments passed to any reachability function
    args = []

    # coefficients matrix
    A = S.s.A
    push!(args, A)

    # determine analysis mode (sparse/dense) for lazy_expm mode
    if A isa SparseMatrixExp
        push!(args, Val(assume_sparse))
    end

    n = statedim(S)
    blocks = kwargs_dict[:blocks]
    partition = convert_partition(kwargs_dict[:partition])

    # Cartesian decomposition of the initial set
    if length(partition) == 1 && length(partition[1]) == n
        info("- No decomposition of X0 needed")
        Xhat0 = [S.x0]
    else
        info("- Decomposing X0")
        tic()
        if lazy_X0
            Xhat0 = S.x0
        elseif !isempty(kwargs_dict[:block_types_init])
            Xhat0 = array(decompose(S.x0, ε=ε_init,
                                    block_types=kwargs_dict[:block_types_init]))
        elseif set_type_init == LazySets.Interval
            Xhat0 = array(decompose(S.x0, set_type=set_type_init, ε=ε_init,
                                    blocks=ones(Int, n)))
        else
            Xhat0 = array(decompose(S.x0, set_type=set_type_init, ε=ε_init))
        end
        tocc()
    end

    # shortcut if only the initial set is required
    if N == 1
        if length(blocks) == 1
            return [Xhat0[blocks[1]]]
        else
            return Xhat0
        end
    end
    push!(args, Xhat0)

    # inputs
    if !assume_homogeneous && inputdim(S) > 0
        U = inputset(S)
    else
        U = nothing
    end
    push!(args, U)

    # overapproximation function for states
    template_directions_symbol = kwargs_dict[:template_directions_iter]
    if template_directions_symbol == :box
        dir = Approximations.BoxDirections
    elseif template_directions_symbol == :oct
        dir = Approximations.OctDirections
    elseif template_directions_symbol == :boxdiag
        dir = Approximations.BoxDiagDirections
    else
        dir = nothing
    end
    if dir != nothing
        overapproximate_fun = (i, x) -> overapproximate(x, dir(length(partition[i])))
    elseif haskey(kwargs_dict, :block_types_iter)
        block_types_iter = block_to_set_map(kwargs_dict[:block_types_iter])
        overapproximate_fun = (i, x) -> (block_types_iter[i] == HPolygon) ?
                                        overapproximate(x, HPolygon, ε_iter) :
                                        overapproximate(x, block_types_iter[i])
    elseif ε_iter < Inf
        overapproximate_fun =
            (i, x) -> overapproximate(x, set_type_iter, ε_iter)
    else
        overapproximate_fun = (i, x) -> overapproximate(x, set_type_iter)
    end
    push!(args, overapproximate_fun)

    # overapproximate function for inputs
    lazy_inputs_interval = kwargs_dict[:lazy_inputs_interval]
    if lazy_inputs_interval == 0
        overapproximate_inputs_fun = (k, i, x) -> overapproximate_fun(i, x)
    else
        # first set
        function _f(k, i, x::LinearMap{MN, NUM}) where {MN, NUM}
            return LazySets.CacheMinkowskiSum(LazySet{NUM}[x])
        end
        # further sets
        function _f(k, i, x::MinkowskiSum)
            @assert (x.X isa CacheMinkowskiSum) "expected CacheMinkowskiSum type, got $(typeof(x.X))"
            push!(array(x.X), x.Y)
            if k % lazy_inputs_interval == 0
                return overapproximate_fun(i, x.X)
            end
            return x.X
        end
        overapproximate_inputs_fun = _f
    end
    push!(args, overapproximate_inputs_fun)

    # ambient dimension
    push!(args, n)

    # number of computed sets
    push!(args, N)

    # preallocate output vector and add mode-specific block(s) argument
    if algorithm == "explicit"
        push!(args, blocks)
        push!(args, partition)
        res = Vector{CartesianProductArray{numeric_type}}(N)
        algorithm_backend = "explicit_blocks"
    else
        error("Unsupported algorithm: ", algorithm)
    end
    push!(args, res)

    # call the adequate function with the given arguments list
    info("- Computing successors")
    tic()
    available_algorithms[algorithm_backend]["func"](args...)
    tocc()

    # return the result
    return res
end
