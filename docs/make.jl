ENV["GKSwstype"] = "100"  # set 'GR environment' to 'no output' (for Travis CI)
using Documenter, Reachability

DocMeta.setdocmeta!(Reachability, :DocTestSetup, :(using Reachability); recursive=true)

makedocs(
    sitename = "Reachability.jl",
    modules = [Reachability],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/juliareach.css"]),
    pages = [
        "Home" => "index.md",
        "Library" => Any[
        "User interface" => "lib/interface.md",
        "Systems" => "lib/systems.md",
        "Algorithms" => "lib/algorithms.md",
        "Transformations" => "lib/transformations.md",
        "Discretization" => "lib/discretize.md",
        "Distributed computations" => "lib/distributed.md"],
        "Publications" => "publications.md",
        "Citations" => "citations.md",
        "About" => "about.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaReach/Reachability.jl.git"
)
