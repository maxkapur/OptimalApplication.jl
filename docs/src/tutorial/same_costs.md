# Working with same-costs markets

The same-costs market is a special case of the college application problem in which all schools have the same application cost. In this case, the budget constraint reduces to a limit on the number of schools to which you can apply:

```math
    \begin{align*}
    \text{maximize}\quad & v(\mathcal{X}) = \sum_{j\in \mathcal{X}} \Bigl( f_j t_j \prod_{\substack{i \in \mathcal{X}: \\ i > j}} (1 - f_{i}) \Bigr) \\
    \text{subject to}\quad & \mathcal{X} \subseteq \mathcal{C}, ~~\lvert \mathcal{X} \rvert \leq h
    \end{align*}
```

Let's consider a random instance of this problem with ``m = 15`` schools and ``h = 5``.

```@example 1
using OptimalApplication
f = rand(15)
t = rand(15)
mkt = SameCostsMarket(f, t, 5)
```

# Enumeration solver

We can solve this problem using [`optimalportfolio_enumerate`](@ref):

```@example 1
X, v = optimalportfolio_enumerate(mkt)
```
This function returns the optimal portfolio as a vector `X` of school indices, and its valuation `v`.

!!! warning

    Enumeration algorithms have exponential time complexity and should only be used in small instances with `m ≤ 20`.
    Their intended use is to test the validity of more-efficient solvers.

# [The nestedness property](@id nestednessproperty)

The optimal solutions for this special case of the college application problem satisfy a property called the
*nestedness property.* When ``\mathcal{X}_h`` denotes the optimal solution for a given market when the limit is ``h``,
we have
```math
\mathcal{X}_1 \subset \mathcal{X}_2  \subset \dots \subset \mathcal{X}_m.
```
This means that there exists a permutation of the schools `X` such that `X[1:h]` gives the optimal portfolio of
size ``h`` for all ``h``.

# Efficient solvers

The functions [`applicationorder_list`](@ref) and [`applicationorder_heap`](@ref) produces this permutation `X`
as well as a vector `V` indicating the utility associated with each optimum. When `mkt.h < mkt.m`, only the first
`mkt.h` entries are computed. To verify the correctness of the nestedness property, however, let's make a copy
of our market with `mkt.h == 15` and check that the first ``h = 5`` entries agree (up to permutation).

```@example 1
mkt_long = SameCostsMarket(f, t, 15)
X_long, V_long = applicationorder_list(mkt)
X_long[1:5], X
```

The valuation should also be the same:

```@example 1
V_long[5], v
```
The output of `applicationorder_heap` is identical to that of `applicationorder_list` but the function uses a different
internal data structure. In the vast majority of cases, the list version is faster.
