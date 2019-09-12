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
        "Home" => Any["Index" => "index.md",
                      "Publications" => "publications.md",
                      "Citations" => "citations.md",
                      "About" => "about.md"],

        "Tutorials" => Any["Linear systems" => "man/linear.md",
                           "Affine systems" => "man/affine.md",
                           "Nonlinear systems" => "man/nonlinear.md",
                           "Hybrid systems" => "man/hybrid.md"],

       "Algorithms" => Any["BFFPSV18" => "algo/BFFPSV18.md"],

        "API" => Any["User interface" => "lib/interface.md",
                     "Systems" => "lib/systems.md",
                     "Algorithms" => "lib/algorithms.md",
                     "Properties" => "lib/properties.md",
                     "Transformations" => "lib/transformations.md",
                     "Discretization" => "lib/discretize.md",
                     "Distributed computations" => "lib/distributed.md"]
    ]
)

deploydocs(
    repo = "github.com/JuliaReach/Reachability.jl.git"
)
