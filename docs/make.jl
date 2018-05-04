ENV["GKSwstype"] = "100"  # set 'GR environment' to 'no output' (for Travis CI)
using Documenter, Reachability

makedocs(
    doctest=false, # use this flag to skip doctests (saves time!)
    modules = [Reachability],
    format = :html,
    assets = ["assets/juliareach.css"],
    sitename = "Reachability.jl",
    pages = [
        "Home" => "index.md",
        "Examples" => Any["SLICOT Benchmarks" => "examples/slicot.md",
                          "Hamiltonian Mechanics" => "examples/hamiltonianmechanics.md"],
        "Library" => Any["User interface" => "lib/interface.md",
                         "Systems" => "lib/systems.md",
                         "Transformations" => "lib/transformations.md",
                         "Discretization" => "lib/discretize.md",
                         "Distributed computations" => "lib/distributed.md"],
        "About" => "about.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaReach/Reachability.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps = nothing,
    make = nothing
)
