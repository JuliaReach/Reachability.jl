import LazySets.Approximations:decompose, overapproximate
export solve

"""
    default_operator(system::InitialValueProblem)

Return the default continous post operator for the initial value problem of a
discrete or continuous system.

### Input

- `system` -- an initial value problem wrapping a mathematical system (continuous or
              discrete) and a set of initial states

### Output

A continuous post operator with default options.
"""
function default_operator(system::InitialValueProblem)
    if islinear(system)
        op = BFFPSV18()
    else
        error("no default reachability algorithm available for system of " *
              "type $(typeof(system))")
    end
    return op
end

"""
    solve(system, options)  or  solve(system, :key1 => val1, [...], keyK => valK)

Solves a reachability problem s.t. the given options.
If some options are not defined, we may fall back to default values.

### Input

- `system`    -- a (discrete or continuoues) system specification
- `options`   -- algorithm options for solving the problem
- `algorithm` -- (optional, default: dispatched on the system's type) the
                 reachability algorithm for the computation

### Output

A solution object whose content depends on the input options.

### Notes

To see all available input options, see
`keys(Reachability.available_keywords.dict)`.
"""
function solve(system::InitialValueProblem,
               options::Options;
               op::ContinuousPost=default_operator(system))
    solve!(system, copy(options), op=op)
end

solve(system::AbstractSystem, options::Pair{Symbol,<:Any}...) =
    solve(system, Options(Dict{Symbol,Any}(options)))

function solve!(system::InitialValueProblem{<:Union{AbstractContinuousSystem,
                                                     AbstractDiscreteSystem}},
                options_input::Options;
                op::ContinuousPost=default_operator(system),
                invariant::Union{LazySet, Nothing}=nothing
               )::AbstractSolution
    system = normalize(system)
    options = init!(op, system, options_input)

    # coordinate transformation
    if options[:coordinate_transformation] != ""
        info("Transformation...")
        (system, transformation_matrix) =
            @timing transform(system, options[:coordinate_transformation])
        options[:transformation_matrix] = transformation_matrix
        invariant = options[:coordinate_transformation] * invariant
    end

    post(op, system, invariant, options)
end

"""
    solve(system::InitialValueProblem{<:HybridSystem},
          options::Options)::AbstractSolution

Interface to reachability algorithms for a hybrid system PWA dynamics.

### Input

- `system`  -- hybrid system
- `options` -- options for solving the problem
"""
function solve(system::InitialValueProblem{<:HybridSystem, <:LazySet},
               options::Options,
               opC::ContinuousPost=BFFPSV18(),
               opD::DiscretePost=LazyDiscretePost())::AbstractSolution
    return solve!(distribute_initial_set(system), copy(options), opC, opD)
end

# if the initial states are distributed, no need to call distribute_initial_set(system)

"""
    constrained_dimensions(HS::HybridSystem)::Dict{Int,Vector{Int}}

Return all coordinates which appear in any guard or invariant constraint for each location.

### Input

- `HS`  -- hybrid system
"""
function constrained_dimensions(HS::HybridSystem)::Dict{Int,Vector{Int}}
    result = Dict{Int,Vector{Int}}()
    sizehint!(result, nstates(HS))
    for mode in states(HS)
        vars = Vector{Int}()
        append!(vars, constrained_dimensions(stateset(HS, mode)))
        for transition in out_transitions(HS, mode)
            append!(vars, constrained_dimensions(stateset(HS, transition)))
        end
        result[mode] = unique(vars)
    end

    return result
end

"""
    constrained_blocks(blocks, nonzero_vars)::Dict{Int,Vector{Int}}

Return all blocks which contain nonzero coordinates.

### Input

- `blocks`  -- blocks structure
"""
function constrained_blocks(blocks, nonzero_vars)::Dict{Int,Vector{Int}}
    result = Dict{Int,Vector{Int}}()
    sizehint!(result, nstates(HS))
    for mode in states(HS)
        vars = Vector{Int}()
        append!(vars, constrained_dimensions(stateset(HS, mode)))
        for transition in out_transitions(HS, mode)
            append!(vars, constrained_dimensions(stateset(HS, transition)))
        end
        result[mode] = unique(vars)
    end

    return result
end


function init_states_sys_from_init_set_sys(
        system::InitialValueProblem{<:HybridSystem, <:LazySet{N}}) where N<:Real
    HS = system.s
    X0 = system.x0
    inits = [(n, X0) for n in states(HS)]
    return InitialValueProblem(HS, inits)
end

function solve(system::InitialValueProblem{<:HybridSystem,
                                           <:Vector{<:Tuple{Int64,<:LazySet}}},
               options::Options,
               opC::ContinuousPost=BFFPSV18(),
               opD::DiscretePost=LazyDiscretePost())::AbstractSolution
    return solve!(system, copy(options), opC, opD)
end

function solve!(system::InitialValueProblem{<:HybridSystem,
                                            <:Vector{<:Tuple{Int64,<:LazySet{N}}}},
               options_input::Options,
               opC::ContinuousPost,
               opD::DiscretePost
              )::AbstractSolution where N<:Real
    # update global variable
    global discrete_post_operator = opD

    HS = system.s
    init_sets = system.x0
    options = init!(opD, HS, options_input)
    time_horizon = options[:T]
    max_jumps = options[:max_jumps]
    inout_map = nothing
    property = options[:mode] == "check" ? options[:property] : nothing
    nonzero_vars = constrained_dimensions(HS)

    # waiting_list entries:
    # - (discrete) location
    # - (set of) continuous-time reach sets
    # - number of previous jumps
    waiting_list = Vector{Tuple{Int, ReachSet{LazySet{N}, N}, Int}}()

    for (loc_id, x0) in init_sets
        loc = HS.modes[loc_id]
        source_invariant = loc.X

        # TODO temporary conversion
        if source_invariant isa HalfSpace
            source_invariant = HPolyhedron([source_invariant])
        end

        #project source invariant and x0
        #TODO add an option for low-dim
        source_invariant_proj = LazySets.Approximations.project(source_invariant, nonzero_vars[loc_id], LinearMap)
        x0_proj = LazySets.Approximations.project(x0, nonzero_vars[loc_id], LinearMap)
        low_dim_x0_inter = overapproximate(source_invariant_proj ∩ x0_proj, Hyperrectangle)
        if !isempty(low_dim_x0_inter)
            loc_x0set = intersection(source_invariant, x0)
            if !isempty(loc_x0set)
                push!(waiting_list, (loc_id,
                    ReachSet{LazySet{N}, N}(loc_x0set, zero(N), zero(N)), 0))
            end
        end
    end

    # passed_list maps the (discrete) location to the (set of) continuous-time
    # reach sets
    passed_list = options[:fixpoint_check] == :none ?
        nothing :
        Vector{Vector{ReachSet{LazySet{N}, N}}}(undef, nstates(HS))

    Rsets = Vector{ReachSet{LazySet{N}, N}}()
    while (!isempty(waiting_list))
        loc_id, X0, jumps = pop!(waiting_list)
        loc = HS.modes[loc_id]
        source_invariant = loc.X

        # compute reach tube
        options_copy = copy(options)
        options_copy_low = copy(options)
        options_copy.dict[:T] = time_horizon - X0.t_start
        options_copy.dict[:project_reachset] = false
        if property != nothing
            # TODO temporary hack, to be resolved in #447
            options_copy[:mode] = "reach"
        end
        if haskey(options_copy, :block_types) &&
                options_copy.dict[:block_types] == nothing
            delete!(options_copy.dict, :block_types)
        end
        if haskey(options_copy, :blocks)
            delete!(options_copy.dict, :blocks)
        end
        #TODO check low-dim option
        options_copy_low = copy(options_copy)
        options_copy_low.dict[:vars] = nonzero_vars[loc_id]
        reach_tube = solve!(IVP(loc, X0.X),
                            options_copy,
                            op=opC,
                            invariant=source_invariant)
        inout_map = reach_tube.options[:inout_map]  # TODO temporary hack
        # get the property for the current location
        property_loc = property isa Dict ?
                       get(property, loc_id, nothing) :
                       property
        if property_loc != nothing
            for (i, reach_set) in enumerate(reach_tube.Xk)
                if !check(property_loc, reach_set.X)
                    return CheckSolution(false, i, options)
                end
            end
        end

        # add the very first initial approximation
        if passed_list != nothing &&
                (!isassigned(passed_list, loc_id) || isempty(passed_list[loc_id]))
            reach_set = reach_tube.Xk[1]
            # TODO For lazy X0 the fixpoint check is likely to fail, so we
            # currently ignore that. In general, we want to add an
            # *underapproximation*, which we currently do not support.
            @assert reach_set.X isa CartesianProductArray
            X0hat = array(reach_set.X)[1]
            if !(X0hat isa ConvexHull)
                Xoa = LazySets.Approximations.overapproximate(reach_set.X)
                ti, tf = reach_set.t_start, reach_set.t_end
                passed_list[loc_id] = [ReachSet{LazySet{N}, N}(Xoa, ti, tf)]
            end
        end
        low_dim_sets = Vector{ReachSet{LazySet{N}, N}}()
        #TODO Add option for lowdim
        if (opD isa LazyLowDimDiscretePost)
            #println("LazyLowDimDiscretePost")

            low_dim_sets = Vector{ReachSet{LazySet{N}, N}}()
        # count_Rsets counts the number of new reach sets added to Rsets

            count_Rsets = tube⋂inv!(opD, reach_tube.Xk, loc.X, Rsets, low_dim_sets,
                                [X0.t_start, X0.t_end], nonzero_vars[loc_id])

            if jumps == max_jumps
                continue
            end

            post(opD, HS, waiting_list, passed_list, loc_id, Rsets, low_dim_sets, count_Rsets,
                 jumps, nonzero_vars[loc_id], options)

        else
        #    println("!LazyLowDimDiscretePost")

            count_Rsets = tube⋂inv!(opD, reach_tube.Xk, loc.X, Rsets,
                                [X0.t_start, X0.t_end])
            if jumps == max_jumps
                continue
            end

            post(opD, HS, waiting_list, passed_list, loc_id, Rsets, count_Rsets,
                 jumps, options)

        end
        #println("Discrete post ends")

    end
    if options[:mode] == "check"
        return CheckSolution(true, -1, options)
    end

    # Projection
    options[:inout_map] = inout_map
    if options[:project_reachset] || options[:projection_matrix] != nothing
        info("Projection...")
        RsetsProj = @timing project(Rsets, options)
    else
        RsetsProj = Rsets
    end
    return ReachSolution(RsetsProj, options)
end
