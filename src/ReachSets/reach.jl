import LazySets.CacheMinkowskiSum

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
function reach(S::AbstractDiscreteSystem,
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
    dir = interpret_template_direction_symbol(
        kwargs_dict[:template_directions_init])
    block_sizes = compute_block_sizes(partition)

    # Cartesian decomposition of the initial set
    if length(partition) == 1 && length(partition[1]) == n
        info("- No decomposition of X0 needed")
        Xhat0 = [S.x0]
    else
        info("- Decomposing X0")
        tic()
        if lazy_X0
            Xhat0 = array(decompose_helper(S.x0, block_sizes, n))
        elseif dir != nothing
            Xhat0 = array(decompose(S.x0, directions=dir, blocks=block_sizes))
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
    dir = interpret_template_direction_symbol(
        kwargs_dict[:template_directions_iter])
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
    if lazy_inputs_interval == nothing
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

    # output function: linear map with the given matrix
    output_function = kwargs_dict[:output_function] != nothing ?
        (x -> kwargs_dict[:output_function] * x) :
        nothing
    push!(args, output_function)

    # preallocate output vector and add mode-specific block(s) argument
    push!(args, blocks)
    push!(args, partition)
    if output_function == nothing
        res = Vector{CartesianProductArray{numeric_type}}(N)
    else
        res = Vector{Hyperrectangle{numeric_type}}(N)
    end

    # choose algorithm backend
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
    tic()
    available_algorithms[algorithm_backend]["func"](args...)
    tocc()

    # return the result
    return res
end

function reach(S::AbstractContinuousSystem,
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

    # ===================
    # Time discretization
    # ===================
    if system isa InitialValueProblem{<:AbstractContinuousSystem}
        info("Time discretization...")
        tic()
        Δ = discretize(
            system,
            options[:δ],
            approx_model=options[:approx_model],
            pade_expm=options[:pade_expm],
            lazy_expm=options[:lazy_expm_discretize],
            lazy_sih=options[:lazy_sih]
            )
        tocc()
    else
        Δ = system
    end
    

    reach(Δ, N, algorithm=algorithm, ε_init=ε_init, set_type_init=set_type_init,
          ε_iter=ε_iter, assume_sparse=assume_sparse, assume_homogeneous=assume_homogeneous,
          numeric_type=numeric_type, lazy_X0=lazy_X0, kwargs...)
end


function matrix_conversion(Δ, options)

    # ===================================
    # Sparse/dense/lazy matrix conversion
    # ===================================
    A = Δ.s.A
    create_new_system = false
    if !options[:lazy_expm] && options[:lazy_expm_discretize]
        # convert SparseMatrixExp to eplicit matrix
        info("Making lazy matrix exponential explicit...")
        tic()
        n = options.dict[:n]
        if options[:assume_sparse]
            B = sparse(Int[], Int[], eltype(A)[], n, n)
        else
            B = Matrix{eltype(A)}(n, n)
        end
        for i in 1:n
            B[i, :] = get_row(A, i)
        end
        A = B
        create_new_system = true
        tocc()
    end
    if options[:assume_sparse]
        if A isa SparseMatrixExp
            # ignore this case
        elseif !method_exists(sparse, Tuple{typeof(A)})
            info("`assume_sparse` option cannot be applied to a matrix of " *
                 "type $(typeof(A)) and will be ignored")
        elseif !(A isa AbstractSparseMatrix)
            # convert to sparse matrix
            A = sparse(A)
            create_new_system = true
        end
    end
    if create_new_system
        # set new matrix
        if method_exists(inputset, Tuple{typeof(Δ.s)})
            Δ = DiscreteSystem(A, Δ.x0, inputset(Δ))
        else
            Δ = DiscreteSystem(A, Δ.x0)
        end
    end
    return Δ
end
 