import LazySets.CacheMinkowskiSum

"""
    check_property(S, options)

Interface to property checking algorithms for an LTI system.

### Input

- `S`                  -- LTI system, discrete or continuous
- `options`            -- additional options

### Notes

A dictionary with available algorithms is available via
`Properties.available_algorithms`.
"""
function check_property(S::IVP{<:AbstractDiscreteSystem},
                        options::Options
                       )::Int
    # list containing the arguments passed to any reachability function
    args = []

    # coefficients matrix
    A = S.s.A
    push!(args, A)

    # determine analysis mode (sparse/dense) for lazy_expm mode
    if A isa SparseMatrixExp
        push!(args, Val(options[:assume_sparse]))
    end

    n = statedim(S)
    blocks = options[:blocks]
    partition = convert_partition(options[:partition])
    dir = interpret_template_direction_symbol(
        options[:template_directions_init])
    block_sizes = compute_block_sizes(partition)
    N = options[:N]
    ε_init = options[:ε_init]
    set_type_init = options[:set_type_init]
    set_type_iter = options[:set_type_iter]

    # Cartesian decomposition of the initial set
    if length(partition) == 1 && length(partition[1]) == n
        info("- No decomposition of X0 needed")
        Xhat0 = [S.x0]
    else
        info("- Decomposing X0")
        tic()
        if options[:lazy_X0]
            Xhat0 = array(decompose_helper(S.x0, block_sizes, n))
        elseif dir != nothing
            Xhat0 = array(decompose(S.x0, directions=dir, blocks=block_sizes))
        elseif !isempty(options[:block_types_init])
            Xhat0 = array(decompose(S.x0, ε=ε_init,
                                    block_types=options[:block_types_init]))
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
            Xhat0_mod = Xhat0[blocks[1]]
        else
            Xhat0_mod = CartesianProductArray(Xhat0)
        end
        return check_property(Xhat0_mod, options[:property]) ? 0 : 1
    end
    push!(args, Xhat0)

    # inputs
    if !options[:assume_homogeneous] && inputdim(S) > 0
        U = inputset(S)
    else
        U = nothing
    end
    push!(args, U)

    # raw overapproximation function
    dir = interpret_template_direction_symbol(
        options[:template_directions_iter])
    if dir != nothing
        overapproximate_fun =
            (i, x) -> overapproximate(x, dir(length(partition[i])))
    elseif haskey(options.dict, :block_types_iter)
        block_types_iter = block_to_set_map(options[:block_types_iter])
        overapproximate_fun = (i, x) -> block_types_iter[i] == HPolygon ?
                              overapproximate(x, HPolygon, set_type_init) :
                              overapproximate(x, block_types_iter[i])
    elseif set_type_init < Inf
        overapproximate_fun =
            (i, x) -> overapproximate(x, set_type_iter, set_type_init)
    else
        overapproximate_fun = (i, x) -> overapproximate(x, set_type_iter)
    end

    # overapproximate function for inputs
    lazy_inputs_interval = options[:lazy_inputs_interval]
    if lazy_inputs_interval == nothing
        overapproximate_inputs_fun = (k, i, x) -> overapproximate_fun(i, x)
    else
        @assert lazy_inputs_interval isa Function "illegal internal value " *
            "$lazy_inputs_interval for option :lazy_inputs_interval"
        # first set in a series
        function _f(k, i, x::LinearMap{NUM}) where {NUM}
            @assert k == 1 "a LinearMap is only expected in the first iteration"
            return CacheMinkowskiSum(LazySet{NUM}[x])
        end
        # further sets of the series
        function _f(k, i, x::MinkowskiSum{NUM, <:CacheMinkowskiSum}) where NUM
            if set_type_init == Inf
                # forget sets if we do not use epsilon-close approximation
                forget_sets!(x.X)
            end
            push!(array(x.X), x.Y)
            if lazy_inputs_interval(k)
                # overapproximate lazy set
                y = overapproximate_fun(i, x.X)
                return CacheMinkowskiSum(LazySet{NUM}[y])
            end
            return x.X
        end
        function _f(k, i, x)
            # other set types
            if lazy_inputs_interval(k)
                # overapproximate lazy set
                return overapproximate_fun(i, x.X)
            end
            return x
        end
        overapproximate_inputs_fun = _f
    end
    push!(args, overapproximate_inputs_fun)

    # ambient dimension
    push!(args, n)

    # number of computed sets
    push!(args, N)

    # add mode-specific block(s) argument
    algorithm = options[:algorithm]
    if algorithm == "explicit"
        push!(args, blocks)
        push!(args, partition)
        algorithm_backend = "explicit_blocks"
    else
        error("Unsupported algorithm: ", algorithm)
    end

    # add eager/lazy checking option
    push!(args, options[:eager_checking])

    # add property
    push!(args, options[:property])

    # call the adequate function with the given arguments list
    info("- Computing successors")
    tic()
    answer = available_algorithms[algorithm_backend]["func"](args...)
    tocc()

    # return the result
    return answer
end

function check_property(S::IVP{<:AbstractContinuousSystem},
                        options::Options
                       )::Int
    # ===================
    # Time discretization
    # ===================
    info("Time discretization...")
    tic()
    Δ = discretize(
        S,
        options[:δ],
        approx_model=options[:approx_model],
        pade_expm=options[:pade_expm],
        lazy_expm=options[:lazy_expm_discretize],
        lazy_sih=options[:lazy_sih]
        )
    tocc()
    check_property(Δ, options)
end
