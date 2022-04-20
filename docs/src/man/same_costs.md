# Same-cost markets

The same-cost market is a special case of the college application problem in which all schools have the same application cost. In this case, the budget constraint reduces to a limit on the number of schools to which you can apply. We maximize
``\operatorname{E}\Bigl[\max\bigr\{0, \max\{t_j Z_j : j \in \mathcal{X}\}\bigr\}\Bigr]``
subject to
``\lvert \mathcal{X} \rvert \leq h .``


## Input and output

```@docs
SameCostsMarket

valuation
```

## Solution algorithms

```@docs
optimalportfolio_enumerate(mkt::SameCostsMarket)
applicationorder_list
applicationorder_heap
```