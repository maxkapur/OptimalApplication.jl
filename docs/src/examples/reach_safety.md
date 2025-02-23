# Reach schools and safety schools

Let's compute the optimal portfolio for a [`VariedCostsMarket`](@ref), plot the results,
and see if we can see any patterns.

First we will read in our packages and generate a random market.

```@example 2
using UnicodePlots, Crayons
using OptimalApplication
mkt = VariedCostsMarket(30)
scatterplot(
    mkt.f,
    mkt.t,
    marker=[k == 10 ? "T" : "$k" for k in mkt.g],
    color=crayon"yellow",
    xlabel="f",
    ylabel="t",
    width=60,
    height=25
)
```

As the plot shows, calling `VariedCostsMarket(m)` generates a random problem in which
`f` and `t` are negatively correlated. This makes the problem more interesting. By default,
the budget `mkt.H` is chosen so that about half the schools end up in the optimal portfolio.
(The numbers in the plot indicate the application cost, where `T` means "ten.")

We can find the optimal solution and objective value using [`optimalportfolio_dynamicprogram`](@ref).

```@example 2
X, v = optimalportfolio_dynamicprogram(mkt)
```

Now let's plot the results:

```@example 2
X_comp = setdiff(1:mkt.m, X)
pl = scatterplot(
    mkt.f[X_comp],
    mkt.t[X_comp],
    marker=[k == 10 ? "T" : "$k" for k in mkt.g[X_comp]],
    color=crayon"red",
    name="don't apply",
    xlabel="f",
    ylabel="t",
    width=60,
    height=25
)
scatterplot!(
    pl,
    mkt.f[X],
    mkt.t[X],
    marker=[k == 10 ? "T" : "$k" for k in mkt.g[X]],
    color=crayon"green",
    name="apply"
)
```

The student should apply to the schools shown in green. For a large budget like this, schools with high utility and
low admissions probability ("reach schools") tend to be preferable to schools with low utility and high
admissions probability ("safety schools") because the number of schools in the portfolio (and therefore the
odds of getting at least one acceptance) is high.

!!! note

    In creating the plot above, `mkt.f[X]` and so on retrieve the data in the correct permutation because
    the random constructor ensures that `mkt.perm = 1:mkt.m`. For a general input data, you need to use
    `mkt.f[invperm(mkt.perm)[X]]`, or simply save `f` and `t` beforehand.
