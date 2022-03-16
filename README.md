# OptimalApplication

Optimal college application strategy with homogeneous and heterogeneous application costs.

Basic usage: 

````julia
julia> using OptimalApplication

julia> mkt = SameCostsMarket(
                # Probability of getting into each school
                [0.2, 0.5, 0.1, 0.6, 0.1],
                # Utility values (must be sorted)
                [1, 4, 5, 7, 8],
                # Number of schools `h` to apply to. By nestedness property, 
                # we can obtain the solution for all `h` by setting `h = m`, 
                # where `m` is the number of schools in the market.
                5
            )
SameCostsMarket(5, [0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 5, 7, 8], 2, [0.2, 2.0, 0.5, 4.2, 0.8], [0.8, 0.5, 0.9, 0.4, 0.9], [1, 2, 3, 4, 5])

julia> x, v = applicationorder(mkt)
([4, 2, 5, 3, 1], [4.2, 5.0, 5.3, 5.4079999999999995, 5.4403999999999995])

julia> x[1:4], v[4] 
([4, 2, 5, 3], 5.4079999999999995)
````

This means that when $h = 4$, the optimal portfolio is $\{4, 2, 5,3\}$ and its valuation is $5.408$. 

Example with varied costs:

````julia
julia> mkt = VariedCostsMarket(
                # Probability of getting into each school
                rand(50),
                # Utility values (must be sorted)
                rand(40:60, 50) |> sort,
                # Cost of applying to each school
                rand(5:10, 50),
                # Budget to spend on applications
                100
            );

julia> x, v = optimalportfolio_dynamicprogram(mkt)
([50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 38, 35, 28], 59.66127736008859)
````

For a large market like this, we may be content with an $\varepsilon$-approximate solution: 

````julia
julia> x, v = optimalportfolio_fptas(mkt, 0.25)
([26, 46, 50], 59.623041055190996)
````

The package also includes `optimalportfolio_enumerate()` and `optimalportfolio_branchbound()`, which are inefficient algorithms of primarily theoretical interest. 
