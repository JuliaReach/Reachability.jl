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
    if isaffine(system)
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

function solve!(problem::InitialValueProblem,
                options::Options;
                op::ContinuousPost=default_operator(problem)
               )

    @assert statedim(problem) == dim(problem.x0) "the state-space dimension should match the " *
    "dimension of the initial states, but they are of size $(statedim(problem)) and $(dim(problem.x0)) respectively"

    # normalize system to a canonical form if needed
    problem = IVP(normalize(problem.s), problem.x0)

    # initialize the algorithm-specific options
    options = init!(op, problem, options)

    # compute a coordinate transformation if needed
    problem, options = transform(problem, options)

    # run the continuous-post operator
    res = post(op, problem, options)

    # undo a coordinate transformation
    res = backtransform(res, options)
end

"""
    solve(system::InitialValueProblem{<:HybridSystem},
          options::Options)

Interface to reachability algorithms for a hybrid system PWA dynamics.

### Input

- `system`  -- hybrid system
- `options` -- options for solving the problem
"""
function solve(system::InitialValueProblem{<:HybridSystem, <:LazySet},
               options::Options,
               opC::ContinuousPost=BFFPSV18(),
               opD::DiscretePost=LazyDiscretePost())
    return solve!(distribute_initial_set(system), copy(options), opC, opD)
end

# if the initial states are distributed, no need to call distribute_initial_set(system)
function solve(system::InitialValueProblem{<:HybridSystem,
                                           <:Vector{<:Tuple{Int64,<:LazySet}}},
               options::Options,
               opC::ContinuousPost=BFFPSV18(),
               opD::DiscretePost=LazyDiscretePost())
    return solve!(system, copy(options), opC, opD)
end

function solve!(system::InitialValueProblem{<:HybridSystem,
                                            <:Vector{<:Tuple{Int64,<:LazySet{N}}}},
               options_input::Options,
               opC::ContinuousPost,
               opD::DiscretePost
              ) where N<:Real
    # update global variable
    global discrete_post_operator = opD

    HS = system.s
    init_sets = system.x0
    options = init!(opD, HS, options_input)
    time_horizon = options[:T]
    max_jumps = options[:max_jumps]
    property = options[:mode] == "check" ? options[:property] : nothing

    # waiting_list entries:
    # - (discrete) location
    # - (set of) continuous-time reach sets
    # - number of previous jumps
    waiting_list = Vector{Tuple{Int, AbstractReachSet{LazySet{N}}, Int}}()

    for (loc_id, x0) in init_sets
        loc = HS.modes[loc_id]
        source_invariant = loc.X

        if x0 ⊆ source_invariant
            loc_x0set = x0
        else
            loc_x0set = intersection(source_invariant, x0)
        end

        if !isempty(loc_x0set)
            aux = ReachSet{LazySet{N}}(loc_x0set, zero(N), zero(N))
            push!(waiting_list, (loc_id, aux, 0))
        end
    end

    # passed_list maps the (discrete) location to the (set of) continuous-time
    # reach sets
    passed_list = options[:fixpoint_check] == :none ?
        nothing :
        Vector{Vector{AbstractReachSet{LazySet{N}}}}(undef, nstates(HS))

    flowpipes = Vector{Flowpipe}()
    while (!isempty(waiting_list))
        loc_id, X0, jumps = pop!(waiting_list)
        loc = HS.modes[loc_id]
        source_invariant = loc.X

        # compute reach tube
        options_copy = copy(options)
        options_copy.dict[:T] = time_horizon - time_start(X0)
        options_copy.dict[:project_reachset] = false
        if property != nothing
            # TODO temporary hack, to be resolved in #447
            options_copy[:mode] = "reach"
        end

        if opC isa BFFPS19
            opC.options.specified[:HS] = HS
            opC.options.specified[:loc_id] = loc_id
            opC.options.specified[:opD] = opD
        end


        reach_tube = solve!(IVP(loc, set(X0)), options_copy, op=opC)
        reach_sets = reach_tube.flowpipes[1].reachsets

        # get the property for the current location
        property_loc = property isa Dict ?
                       get(property, loc_id, nothing) :
                       property

        # add the very first initial approximation
        if passed_list != nothing &&
                (!isassigned(passed_list, loc_id) || isempty(passed_list[loc_id]))
            reach_set = reach_sets[1]
            # TODO For lazy X0 the fixpoint check is likely to fail, so we
            # currently ignore that. In general, we want to add an
            # *underapproximation*, which we currently do not support.
            X = set(reach_set)
            if !(X isa CartesianProductArray) || !(array(X)[1] isa CH)
                Xoa = overapproximate(X)
                ti, tf = time_start(reach_set), time_end(reach_set)
                passed_list[loc_id] = [ReachSet{LazySet{N}}(Xoa, ti, tf)]
            end
        end

        Rsets = tube⋂inv(opD, reach_sets, loc.X, [time_start(X0), time_end(X0)])

        if property_loc != nothing
            if opD isa DecomposedDiscretePost
                temp_vars = opD.options[:temp_vars]
                n_lowdim = length(temp_vars)
                n = dim(set(X0))
                property_loc_lowdim = project(property_loc, temp_vars)
                for (i, reach_set) in enumerate(reach_sets)
                    X = set(reach_set)
                    if (dim(X) == n_lowdim && n_lowdim < n)
                        if !property_loc_lowdim(X)
                            return CheckSolution(false, i, options)
                        end
                    elseif !property_loc(X)
                        return CheckSolution(false, i, options)
                    end
                end
            else
                for (i, reach_set) in enumerate(Rsets)
                    if !property_loc(set(reach_set))
                        return CheckSolution(false, i, options)
                    end
                end
            end
        end

        if jumps == max_jumps
            continue
        end
        post(opD, HS, waiting_list, passed_list, loc_id, Rsets, jumps, options)
        # create vector with concrete set type (needed by ReachSolution)
        push!(flowpipes, Flowpipe([rs for rs in Rsets]))
    end
    if options[:mode] == "check"
        return CheckSolution(true, -1, options)
    end

    # Projection
    if options[:project_reachset] || options[:projection_matrix] != nothing
        info("Projection...")
        RsetsProj = @timing project(ReachSolution(flowpipes, options))
    else
        RsetsProj = ReachSolution(flowpipes, options)
    end
    return RsetsProj
end
