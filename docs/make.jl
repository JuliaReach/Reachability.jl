using Documenter, Reachability

makedocs(
    modules = [Reachability],
    format = :html,
    sitename = "Reachability.jl",
    pages = [
        "Home" => "index.md",
        "Library" => Any[
        "User interface" => "lib/interface.md"],
        "About" => "about.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaReach/Reachability.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    make = nothing
)