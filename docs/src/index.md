# OptimalApplication

... is a Julia package for solving the *college application problem,* which is the following optimization problem:

```math
    \begin{align*}
    \text{maximize}\quad & v(\mathcal{X}) =  \operatorname{E}\Bigl[\max\bigr\{0,
    \max\{t_j Z_j : j \in \mathcal{X}\}\bigr\}\Bigr] \\
    \text{subject to}\quad & \mathcal{X} \subseteq \mathcal{C}, ~~\sum_{j\in \mathcal{X}} g_j \leq H
    \end{align*}
```

Here ``\mathcal{C} = \{ 1 \dots m\}`` is a collection of schools, ``t_j`` is the utility associated with
attending school ``j``, and ``Z_j`` is a random Bernoulli variable with probability ``f_j``. The problem is
to decide on a subset of schools ``\mathcal{X}`` to apply to, given ``g_j``, the application fee for school
``j``, and ``H`` the total budget to spend on applications. ``v(\mathcal{X})`` represents the expected utility
associated with the application portfolio ``\mathcal{X}``, and by indexing the schools such that each
``t_j \leq t_{j+1}``, it can be written in the following closed form: 

```math
    v(\mathcal{X}) = \sum_{j\in \mathcal{X}} \Bigl( f_j t_j \prod_{\substack{i \in \mathcal{X}: \\ i > j}} (1 - f_{i}) \Bigr)
```

The general college application problem is referred to in this package as the *varied-costs* college application problem.
A special case, called the *same-costs* problem, occurs when each ``g_j = 1``. Then the budget constraint becomes a simple
cardinality constraint ``\lvert \mathcal{X} \rvert \leq h``, where ``h`` is written in lowercase to reflect the change in 
meaning. 

The same-costs college application problem is solvable in polynomial time by a greedy algorithm, whereas the varied-costs
problem is NP-complete (as can be shown by a reduction from the knapsack problem).

This package provides:

- Types [`SameCostsMarket`](@ref) and [`VariedCostsMarket`](@ref) for holding problem instance data.
- Various exact algorithms for `SameCostsMarket`s
- Exact algorithms for `VariedCostsMarket`s
- An approximation algorithm for `VariedCostsMarket`s
- Heuristic algorithms for `VariedCostsMarket`s

```@contents
Pages = [
    "man/same_costs.md",
    "man/varied_costs.md",
]
Depth = 1
```
