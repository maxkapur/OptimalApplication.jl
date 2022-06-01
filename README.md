# OptimalApplication.jl

[![In Development](https://img.shields.io/badge/docs-dev-blue.svg)](https://misc.maxkapur.com/OptimalApplication.jl/dev/)

Optimal college application strategy with homogeneous and heterogeneous application costs.

Basic usage: 

````julia
julia> using OptimalApplication

julia> mkt = SameCostsMarket(
                # Probability of getting into each school
                [0.2, 0.5, 0.1, 0.6, 0.1],
                # Utility values
                [1, 4, 9, 1, 8],
                # Number of schools `h` to apply to. By nestedness property, 
                # we can obtain the solution for all `h` by setting `h = m`, 
                # where `m` is the number of schools in the market.
                5
            );

julia> x, v = applicationorder_list(mkt)
([2, 3, 5, 4, 1], [2.0, 2.7, 3.24, 3.483, 3.5154])

julia> x[1:4], v[4] 
([2, 3, 5, 4], 3.483)
````

This means that when `h = 4`, the optimal portfolio is `{2, 3, 5, 4}` and its valuation is `5.408`.

The function `applicationorder_heap()` works the same way but uses a different internal data structure and may be faster for certain instances. 

Example with varied costs:

````julia
julia> mkt = VariedCostsMarket(
                # Probability of getting into each school
                rand(50),
                # Utility values
                rand(40:60, 50),
                # Cost of applying to each school
                rand(5:10, 50),
                # Budget to spend on applications
                100
            );

julia> x, v = optimalportfolio_dynamicprogram(mkt)
([50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 38, 35, 28], 59.66127736008859)
````

For a large market like this, we may be content with an ε-approximate solution: 

````julia
julia> x, v = optimalportfolio_fptas(mkt, 0.25)
([26, 46, 50], 59.623041055190996)
````

The value `ε = 0.25` means that the value is guaranteed to be no worse than `1 - 0.25` times the value of the optimal portfolio. As we can see above, the typical optimality gap is often much tighter.

The package also includes `optimalportfolio_enumerate()` and `optimalportfolio_branchbound()`, which are inefficient algorithms of primarily theoretical interest.

Finally, we can check the expected value of an arbitrary portfolio using `valuation()`:

````julia
julia> valuation([20, 16, 35], mkt)
35.33785503577366
````

## arXiv paper

If you found this package useful, please consider citing [our arXiv paper](https://arxiv.org/abs/2205.01869). Rolling updates to the paper can be found on the companion repository [CollegeApplication](https://github.com/maxkapur/CollegeApplication).
