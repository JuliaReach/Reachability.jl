"""
    tune_δ(algorithm, time_horizon, [precision], [δ], [options])

Find the threshold step size for a given property.

### Input

- `algorithm`    -- the function that is used to evaluate ``δ``;
                    typically this is the compute function from a model script
- `time_horizon` -- time horizon
- `precision`    -- precision of result (s. a.);
                      default: 1e-5
- `δ`            -- (optional) initial guess for the time step;
                      default: 1e-3
- `options`      -- (optional) options passed to ``algorithm``;
                      default: hylaaOptions

### Algorithm

Finds a value for ``δ`` such that `algorithm`:

  1. Returns `true` for ``δ`` and
  2. Returns `false` for (``δ`` + ``precision``).
"""
function tune_δ(algorithm::Function,
                time_horizon::Float64,
                precision::Float64=1e-5,
                δ::Float64=1e-3,
                options::Options=hylaaOptions
               )::Float64
    biggest_working = 0.
    smallest_failing = time_horizon
    smallest_failing_modified = false

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

        assert((δ > biggest_working && δ < smallest_failing) || (!smallest_failing_modified && δ == smallest_failing))

        # choose N according to δ
        N = ceil(Int64, time_horizon / δ)

        # check property
        tic()
        answer = algorithm(N, δ, options)
        runtime = toq()
        if answer
            # success
            biggest_working = δ
            runtime_best = runtime

            if δ >= time_horizon
                # exceeded time horizon -> terminate
                break
            end

            # refine δ
            if !smallest_failing_modified
                # no upper bound yet, choose δ that is 50 % bigger
                δ_new = min(δ * 1.5, time_horizon)
            else
                if step <= 0.
                    # previous iteration was no increase, choose new step
                    # choose step that 10 % closer in the uncertainty interval
                    step = (smallest_failing - δ) * 0.1
                end
                δ_new = δ + step
                while δ_new >= smallest_failing
                    # reached upper bound, reduce step by 50 %
                    step = step * 0.5
                    δ_new = δ + step
                end
                assert(step > 0.)
            end
            assert(δ_new > δ)
        else
            # failure
            smallest_failing = δ
            smallest_failing_modified = true

            # refine δ
            if biggest_working == 0.
                # no lower bound yet, choose δ that is 50 % smaller
                δ_new = δ * 0.5
            else
                if step >= 0.
                    # previous iteration was no decrease, choose new step
                    # choose step that 10 % closer in the uncertainty interval
                    step = (biggest_working - δ) * 0.1
                end
                δ_new = δ + step
                while δ_new <= biggest_working
                    # reached upper bound, reduce step by 50 %
                    step = step * 0.5
                    δ_new = δ + step
                end
                assert(step < 0.)
            end
            assert(δ_new < δ)
        end
        if smallest_failing - biggest_working <= precision
            # reached precision window -> terminate
            break
        end

        δ = δ_new
    end
    info("biggest working δ found: ", biggest_working, " (runs in ", runtime_best, " sec)")
    info("smallest failing δ found: ", smallest_failing_modified ? smallest_failing : "--")
    return δ
end
