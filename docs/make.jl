using Documenter
using OptimalApplication

DocMeta.setdocmeta!(OptimalApplication, :DocTestSetup, :(using OptimalApplication); recursive=true)

makedocs(
    sitename="OptimalApplication",
    modules=[OptimalApplication],
    format=Documenter.HTML(ansicolor=true),
    # format = Documenter.HTML(prettyurls = false), # For local browsing
    pages=Any[
        "Home"=>"index.md",
        "Tutorial"=>[
            "Market I/O" => "tutorial/market_io.md",
            "Same-costs problem" => "tutorial/same_costs.md",
            "Varied-costs problem" => "tutorial/varied_costs.md",
        ],
        "Public API"=>[
            "Market I/O" => "public_api/market_io.md",
            "Enumeration solvers" => "public_api/enumeration_solvers.md",
            "Fast, exact solvers" => "public_api/fast_solvers.md",
            "Approximation methods" => "public_api/approximation.md",
            "Heuristic solvers" => "public_api/heuristics.md",
        ],
        "Examples"=>[
            "Valuation curve" => "examples/valuation.md",
            "Reach & safety schools" => "examples/reach_safety.md",
        ],
    ],
)
