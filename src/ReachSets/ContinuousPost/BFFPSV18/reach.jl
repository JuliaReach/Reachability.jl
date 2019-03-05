import LazySets: CacheMinkowskiSum,
                 isdisjoint

import ..Utils: LDS, CLCDS

"""
    reach(S, invariant, options)

Interface to reachability algorithms for an LTI system.

### Input

- `S`         -- LTI system, discrete or continuous
- `invariant` -- invariant
- `options`   -- additional options

### Output

A sequence of [`ReachSet`](@ref)s.

### Notes

A dictionary with available algorithms is available via
`Reachability.available_algorithms`.

The numeric type of the system's coefficients and the set of initial states
is inferred from the first parameter of the system (resp. lazy set), ie.
`NUM = first(typeof(S.s).parameters)`.
"""
function reach(S::Union{IVP{<:LDS{NUM}, <:LazySet{NUM}},
                        IVP{<:CLCDS{NUM}, <:LazySet{NUM}}},
               invariant::Union{LazySet, Nothing},
               options::TwoLayerOptions
              )::Vector{<:ReachSet} where {NUM <: Real}

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
    N = ceil(Int, options[:T] / options[:δ])
    ε_init = options[:ε_init]
    set_type_init = options[:set_type_init]
    ε_iter = options[:ε_iter]
    set_type_iter = options[:set_type_iter]


    # Cartesian decomposition of the initial set
    if length(partition) == 1 && length(partition[1]) == n
        info("- No decomposition of X0 needed")
        Xhat0 = LazySet{NUM}[S.x0]
    else
        info("- Decomposing X0")
        @timing begin
            if options[:lazy_X0]
                Xhat0 = array(decompose_helper(S.x0, block_sizes, n))
            elseif dir != nothing
                Xhat0 = array(decompose(S.x0, directions=dir,
                                        blocks=block_sizes))
            elseif options[:block_types_init] != nothing &&
                    !isempty(options[:block_types_init])
                Xhat0 = array(decompose(S.x0, ε=ε_init,
                                        block_types=options[:block_types_init]))
            elseif set_type_init == LazySets.Interval
                Xhat0 = array(decompose(S.x0, set_type=set_type_init, ε=ε_init,
                                        blocks=ones(Int, n)))
            else
                Xhat0 = array(decompose(S.x0, set_type=set_type_init, ε=ε_init,
                                        blocks=block_sizes))
            end
        end
    end

    # determine output function: linear map with the given matrix
    output_function = options[:output_function] != nothing ?
        (x -> options[:output_function] * x) :
        nothing

    # preallocate output vector
    if output_function == nothing
        res_type = ReachSet{CartesianProductArray{NUM, LazySet{NUM}}, NUM}
    else
        res_type = ReachSet{Hyperrectangle{NUM}, NUM}
    end
    res = Vector{res_type}(undef, N)

    # shortcut if only the initial set is required
    if N == 1
        res[1] = res_type(
            CartesianProductArray{NUM, LazySet{NUM}}(Xhat0[blocks]),
            zero(NUM), options[:δ])
        return res
    end
    push!(args, Xhat0)

    # inputs
    if !options[:assume_homogeneous] && inputdim(S) > 0
        U = inputset(S)
    else
        U = nothing
    end
    push!(args, U)

    # overapproximation function for states
    dir = interpret_template_direction_symbol(
        options[:template_directions_iter])
    if dir != nothing
        overapproximate_fun = (i, x) -> overapproximate(x, dir(length(partition[i])))
    elseif options[:block_types_iter] != nothing
        block_types_iter = block_to_set_map(options[:block_types_iter])
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
    lazy_inputs_interval = options[:lazy_inputs_interval]
    if lazy_inputs_interval == lazy_inputs_interval_always
        overapproximate_inputs_fun = (k, i, x) -> overapproximate_fun(i, x)
    else
        # first set in a series
        function _f(k, i, x::LinearMap{NUM}) where {NUM}
            @assert k == 1 "a LinearMap is only expected in the first iteration"
            return CacheMinkowskiSum(LazySet{NUM}[x])
        end
        # further sets of the series
        function _f(k, i, x::MinkowskiSum{NUM, <:CacheMinkowskiSum}) where NUM
            if ε_iter == Inf
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

    # output function
    push!(args, output_function)

    # add mode-specific block(s) argument
    push!(args, blocks)
    push!(args, partition)

    # time step
    push!(args, options[:δ])

    # termination function
    if invariant == nothing
        termination = (k, set, t0) -> termination_N(N, k, set, t0)
    else
        termination =
            (k, set, t0) -> termination_inv_N(N, invariant, k, set, t0)
    end
    push!(args, termination)

    # choose algorithm backend
    algorithm = options[:algorithm]
    if algorithm == "explicit"
        algorithm_backend = "explicit_blocks"
    elseif algorithm == "wrap"
        algorithm_backend = "wrap"
    else
        error("Unsupported algorithm: ", algorithm)
    end
    push!(args, res)

    # call the adequate function with the given arguments list
    info("- Computing successors")
    @timing begin
        index, skip = available_algorithms[algorithm_backend]["func"](args...)
        if index < N || skip
            # shrink result array
            info("terminated prematurely, only computed $index/$N steps")
            deleteat!(res, (skip ? index : index + 1):N)
        end
    end

    # return the result
    return res
end

function reach(system::IVP{<:AbstractContinuousSystem},
               invariant::Union{LazySet, Nothing},
               options::TwoLayerOptions
              )::Vector{<:ReachSet}
    # ===================
    # Time discretization
    # ===================
    info("Time discretization...")
    Δ = @timing discretize(system, options[:δ], algorithm=options[:discretization],
                                                exp_method=options[:exp_method],
                                                sih_method=options[:sih_method])

    Δ = matrix_conversion(Δ, options)
    return reach(Δ, invariant, options)
end

function termination_N(N, k, set, t0)
    return (k >= N, false)
end

function termination_inv_N(N, inv, k, set, t0)
    if k >= N
        return (true, false)
    elseif isdisjoint(set, inv)
        return (true, true)
    else
        return (false, false)
    end
end
