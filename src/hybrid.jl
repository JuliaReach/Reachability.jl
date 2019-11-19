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
               opC::AbstractContinuousPost=BFFPSV18(),
               opD::AbstractDiscretePost=LazyDiscretePost())
    return solve!(distribute_initial_set(system), copy(options), opC, opD)
end

# if the initial states are distributed, no need to call distribute_initial_set(system)
function solve(system::InitialValueProblem{<:HybridSystem,
                                           <:Vector{<:Tuple{Int64,<:LazySet}}},
               options::Options,
               opC::AbstractContinuousPost=BFFPSV18(),
               opD::AbstractDiscretePost=LazyDiscretePost())
    return solve!(system, copy(options), opC, opD)
end

function solve!(system::InitialValueProblem{<:HybridSystem,
                                            <:Vector{<:Tuple{Int64,<:LazySet{N}}}},
               options_input::Options,
               opC::AbstractContinuousPost,
               opD::AbstractDiscretePost
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

    Rsets = Vector{AbstractReachSet{<:LazySet{N}}}()
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

        # get the property for the current location
        property_loc = property isa Dict ?
                       get(property, loc_id, nothing) :
                       property

        # add the very first initial approximation
        if passed_list != nothing &&
                (!isassigned(passed_list, loc_id) || isempty(passed_list[loc_id]))
            reach_set = reach_tube.Xk[1]
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

        # count_Rsets counts the number of new reach sets added to Rsets
        count_Rsets = tube⋂inv!(opD, reach_tube.Xk, loc.X, Rsets,
                                [time_start(X0), time_end(X0)])

        if property_loc != nothing
            if opD isa DecomposedDiscretePost
                temp_vars = opD.options[:temp_vars]
                n_lowdim = length(temp_vars)
                n = dim(set(X0))
                property_loc_lowdim = project(property_loc, temp_vars)
                for (i, reach_set) in enumerate(reach_tube.Xk)
                    X = set(reach_set)
                    if (dim(X) == n_lowdim && n_lowdim < n)
                        if !check(property_loc_lowdim, X)
                            return CheckSolution(false, i, options)
                        end
                    elseif !check(property_loc, X)
                        return CheckSolution(false, i, options)
                    end
                end
            else
                for (i, reach_set) in enumerate(Rsets[length(Rsets) - count_Rsets + 1 : end])
                    if !check(property_loc, set(reach_set))
                        return CheckSolution(false, i, options)
                    end
                end
            end
        end

        if jumps == max_jumps
            continue
        end
        post(opD, HS, waiting_list, passed_list, loc_id, Rsets, count_Rsets,
            jumps, options)

    end
    if options[:mode] == "check"
        return CheckSolution(true, -1, options)
    end

    # create vector with concrete set type (needed by ReachSolution)
    Rsets = [rs for rs in Rsets]

    # Projection
    if options[:project_reachset] || options[:projection_matrix] != nothing
        info("Projection...")
        RsetsProj = @timing project(Rsets, options)
    else
        RsetsProj = Rsets
    end
    return ReachSolution(RsetsProj, options)
end
