# Working with varied-costs markets

The varied-costs market is the general case of the college application problem in which schools
may have different application costs. The optimization problem is as follows:

```math
    \begin{align*}
    \text{maximize}\quad & v(\mathcal{X}) = \sum_{j\in \mathcal{X}} \Bigl( f_j t_j \prod_{\substack{i \in \mathcal{X}: \\ i > j}} (1 - f_{i}) \Bigr) \\
    \text{subject to}\quad & \mathcal{X} \subseteq \mathcal{C}, ~~\sum_{j\in \mathcal{X}} g_j \leq H
    \end{align*}
```

Let's consider a random instance of this problem with ``m = 15`` schools.

```@example 1
using OptimalApplication
f = rand(15)
t = rand(60:80, 15)     # must be integers
g = rand(5:10, 15)      # must be integers
H = sum(g) ÷ 2
mkt = VariedCostsMarket(f, t, g, H)
```

# Enumeration solver

We can solve this problem using [`optimalportfolio_enumerate`](@ref):

```@example 1
optimalportfolio_enumerate(mkt)
```
This function returns the optimal portfolio as a vector `X` of school indices, and its valuation `v`.

A somewhat better kind of enumeration the branch-and-bound scheme, implemented as [`optimalportfolio_branchbound`](@ref):

```@example 1
optimalportfolio_branchbound(mkt)
```

We can see that the results are the same (up to permutation).

!!! warning

    Enumeration algorithms have exponential time complexity and should only be used in small instances with `m ≤ 20`.
    Their intended use is to test the validity of more-efficient solvers.

# Fast solver

OptimalApplication implements a dynamic program that solves this problem exactly in pseudopolynomial-time.
The computation time of [`optimalportfolio_dynamicprogram`](@ref) depends on the values of ``g_j``; smaller is better.

```@example 1
optimalportfolio_dynamicprogram(mkt)
```

# Approximation algorithms

The college application problem admits a fully polynomial-time approximation scheme (FPTAS). Given a market and tolerance
``\varepsilon``, [`optimalportfolio_fptas(mkt, ε)`](@ref): produces a solution that is guaranteed to have an objective value
of at least ``1 - \varepsilon`` times the optimum, in time ``O(m^3 / \varepsilon)``.

```@example 1
optimalportfolio_fptas(mkt, 0.2)
```

!!! warning

    Setting the tolerance `ε` to an extremely small value requires lots of system memory.

# Heuristic algorithms

OptimalApplication also includes two heuristic algorithms. Although they offer no accuracy guarantee, their computation time
is low and they can be quite effective in practice.

[`optimalportfolio_greedy`](@ref) adds schools in descending order by `f_j t_j / g_j` until the budget is exhausted.

```@example 1
optimalportfolio_greedy(mkt)
```

[`optimalportfolio_simulatedannealing`](@ref) starts with the greedy solution, then searches locally according to a randomized
scheme for an improvement.

```@example 1
optimalportfolio_simulatedannealing(mkt)
```
