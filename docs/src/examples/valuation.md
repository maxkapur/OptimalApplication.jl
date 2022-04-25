# Planets example

Let's consider a same-costs market and plot the *valuation curve,* which gives ``v(\mathcal{X}_h)``
as a function of ``h``, where ``\mathcal{X}_h`` is the optimal portfolio for each ``h`` between ``0``
and ``m``.

The market data is as follows:

```math
\begin{align}
f &= (0.39, 0.33, 0.24, 0.24, 0.05, 0.03, 0.1, 0.12) \\
t &= (200, 250, 300, 350, 400, 450, 500, 550)
\end{align}
```

Let's import the `PrettyTables` and `UnicodePlots` libraries for visualizations and read in the market data.


```@example 1
using PrettyTables
using UnicodePlots
using OptimalApplication
f = [0.39, 0.33, 0.24, 0.24, 0.05, 0.03, 0.1, 0.12]
t = [200, 250, 300, 350, 400, 450, 500, 550]
mkt = Market(f, t, 8)
```


Notice how we set ``h = m = 8`` because we know from the [nestedness property](@ref nestednessproperty) that
by calling [`applicationorder_list(mkt)`](@ref), we can get a permutation of the schools that encodes
all the optimal portfolios. Let's do that now:

```@example 1
X, V = applicationorder_list(mkt, verbose=true)
```

``\mathcal{X}_h`` is given by the first ``h`` entries of `X`, and ``v(\mathcal{X}_h)`` is `V[h]`.
Let's plot `V`.

```@example 1
lineplot(0:mkt.m, vcat(0, V), xlabel="h", ylabel="v")
```

This graph is called a *valuation curve.* It gives the applicant's utility as a function
of her budget ``h``. The curve's concave shape is guaranteed by the nestedness property. 



An alternative way of looking at the optimal portfolios is to take the *inverse-permutation* of `X`.
The value of `invperm(X)[j]` gives the *smallest* value of ``h`` for which school ``j`` is in
``\mathcal{X}_h``. This is called the schools *priority* number; if a school has priority number 1,
it is the first school you should apply to, and so on. 

This format lets us create the following table:

```@example 1
priority = invperm(X)
pretty_table(
    Any[f t priority V[priority]], 
    header = ["f", "t", "priority", "valuation"],
    header_crayon = crayon"bold yellow",
)
```


