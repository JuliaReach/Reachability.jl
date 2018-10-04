"""
    tune_δ(algorithm, time_horizon, [precision], [δ], [options])

Find the threshold step size for a given property.

### Input

- `algorithm`        -- the function that is used to evaluate ``δ``;
                        the call should be: `algorithm(N, δ)`
- `time_horizon`     -- time horizon
- `precision`        -- (optional, default: 1e-5) precision of result (s. a.)
- `δ`                -- (optional, default: 1e-1) initial guess for the time
                        step
- `rel_step_size`    -- (optional, default: 0.1) relative step size
- `biggest_working`  -- (optional, default: 0.) biggest known working ``δ``
- `smallest_failing` -- (optional, default: `time_horizon`) smallest known
                        failing ``δ``

### Algorithm

Finds a value for ``δ`` such that `algorithm`

  1. returns `true` for ``δ`` and
  2. returns `false` for (``δ`` + `precision`).
"""
function tune_δ(algorithm::Function,
                time_horizon::Float64,
                precision::Float64=1e-5,
                δ::Float64=1e-1,
                rel_step_size::Float64=0.1;
                biggest_working=0.,
                smallest_failing=time_horizon
                )::Float64
    biggest_working_modified = biggest_working != 0.
    smallest_failing_modified = smallest_failing != time_horizon

    step = 0.
    runtime_best = -1.

    if δ <= 0. || precision <= 0.
        throw(DomainError())
    end

    i = 0
    while true
        i += 1
        info("Iteration ", i)
        info("δ: ", δ)

        @assert (δ >= biggest_working && δ <= smallest_failing) "invalid δ: $δ (not in [$biggest_working, $smallest_failing])"

        # choose N according to δ
        N = ceil(Int, time_horizon / δ)

        # check property
        tic()
        answer = algorithm(N, δ)
        runtime = toq()
        if answer
            # success
            biggest_working = δ
            biggest_working_modified = true
            runtime_best = runtime
        else
            # failure
            smallest_failing = δ
            smallest_failing_modified = true
        end

        # choose new δ
        if smallest_failing - biggest_working <= precision
            # reached precision window -> terminate
            break
        elseif biggest_working_modified && smallest_failing_modified
            # choose δ in the center of the uncertainty interval
            δ_new = smallest_failing + (biggest_working - smallest_failing) / 2.
        elseif biggest_working_modified
            # no upper bound yet, choose δ that is bigger wrt. the relative ball
            δ_new = min(δ * (1. + rel_step_size), time_horizon)
        else
            @assert smallest_failing_modified
            # no lower bound yet, choose δ that is bigger wrt. the relative ball
            δ_new = max(δ * (1. - rel_step_size), 0.)
        end

        # round δ down wrt. precision
        inverse_precision = round(Int, 1. / precision)
        δ_new = floor(δ_new * inverse_precision) / inverse_precision
        if δ_new >= smallest_failing || δ_new <= biggest_working
            # reached precision window -> terminate
            break
        end

        δ = δ_new
    end
    info("biggest working δ found: ", biggest_working, " (runs in ", runtime_best, " sec)")
    info("smallest failing δ found: ", smallest_failing_modified ? smallest_failing : "--")
    return δ
end
