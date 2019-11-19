"""
    default_operator(problem::InitialValueProblem)

Return the default continous post operator for the initial value problem of a
discrete or continuous system.

### Input

- `problem` -- an initial value problem represented by a mathematical system
               and a set of initial states

### Output

A continuous post operator with default options.
"""
function default_operator(problem::InitialValueProblem{})
    if isaffine(problem) || islinear(problem)
        op = BFFPSV18()
    else
        op = TMJets()
    end
    return op
end
