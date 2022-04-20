using Documenter
using OptimalApplication

DocMeta.setdocmeta!(OptimalApplication, :DocTestSetup, :(using OptimalApplication); recursive=true)

makedocs(
    sitename="OptimalApplication",
    modules = [OptimalApplication],
    format = Documenter.HTML(prettyurls = false), # For local browsing
    pages = Any[
            "Home" => "index.md",
            "Manual" => [
                "Same-cost markets" => "man/same_costs.md",
                "Varied-cost markets" => "man/varied_costs.md",
            ],
        ],
)
