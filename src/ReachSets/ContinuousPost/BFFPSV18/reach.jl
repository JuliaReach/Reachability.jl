using ..Utils: LDS, CLCDS

using LazySets: CacheMinkowskiSum,
                 isdisjoint

import LazySets.Approximations: overapproximate

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
    steps = haskey(options.specified, :steps) ? options[:steps] : Vector{Int}()
    n = statedim(S)
    blocks = options[:blocks]
    partition = convert_partition(options[:partition])
    N = ceil(Int, options[:T] / options[:δ])


    # Cartesian decomposition of the initial set
    if length(partition) == 1 && length(partition[1]) == n &&
            options[:block_options_init] == LinearMap
        info("- Skipping decomposition of X0")
        Xhat0 = LazySet{NUM}[S.x0]
    else
        info("- Decomposing X0")
        @timing begin
            Xhat0 = array(decompose(S.x0, options[:partition],
                                    options[:block_options_init]))
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
    if isempty(steps)
       res = Vector{res_type}(undef, N)
    else
       res = Vector{res_type}(undef, length(steps))
    end

    # shortcut if only the initial set is required
    if N == 1
        res[1] = res_type(
            CartesianProductArray{NUM, LazySet{NUM}}(Xhat0[blocks]),
            zero(NUM), options[:δ], 0)
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
    block_options_iter = options[:block_options_iter]
    if block_options_iter isa AbstractVector ||
            block_options_iter isa Dict{Int, Any}
        # individual overapproximation options per block
        overapproximate_fun = (i, X) -> overapproximate(X, block_options_iter[i])
    else
        # uniform overapproximation options for each block
        overapproximate_fun = (i, X) -> overapproximate(X, block_options_iter)
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
            if has_constant_directions(block_options_iter, i)
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
    constraints_list = options[:constraints]
    if invariant == nothing
        termination = (k, set, t0) -> termination_N(N, k, set, t0)
    elseif isempty(constraints_list)
        termination =
            (k, set, t0) -> termination_inv_N(N, invariant, k, set, t0)
    else
        termination =
            (inv, k, set, t0) -> termination_inv_N(N, inv, k, set, t0)
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
    if isempty(constraints_list)
        push!(args, steps)
    else
        push!(args, constraints_list)
    end
    push!(args, invariant)
    # comp_vars = options[:comp_vars]
    # push!(args, comp_vars)

    # call the adequate function with the given arguments list
    info("- Computing successors")
    @timing begin
        index, skip = available_algorithms[algorithm_backend]["func"](args...)
        if isempty(steps) && (index < N || skip)
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
    # println("Time discretization")
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
        return (true, false, inv)
    else
        # println("inv ", inv)
        inv_int = intersection(set, inv)
        if isempty(inv_int)
            return (true, true, inv)
        else
            return (false, false, inv_int)
        end
    end
end

function overapproximate(X::LazySet, pair::Pair)
    return overapproximate(X, pair[1], pair[2])
end

function has_constant_directions(block_options::AbstractVector, i::Int)
    return has_constant_directions(block_options[i], i)
end

function has_constant_directions(block_options::Dict{<:UnionAll, <:Real},
                                 i::Int)
    return has_constant_directions(block_options[i], i)
end

function has_constant_directions(block_options::Pair{<:UnionAll, <:Real},
                                 i::Int)
    return has_constant_directions(block_options[2], i)
end

function has_constant_directions(block_options::Real, i::Int)
    return ε == Inf
end

function has_constant_directions(block_options, i::Int)
    return true
end
