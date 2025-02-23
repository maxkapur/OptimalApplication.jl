# Market input and output

To read in an instance of the same-costs problem, use [`SameCostsMarket(f, t, h)`](@ref):

```@example 1
using OptimalApplication
mkt = SameCostsMarket(
    [0.5, 0.1, 0.9, 0.7],   # admissions probabilities
    [12, 20, 1, 3],         # utility values
    2                       # limit on number of schools to apply to
)
```

The output is an instance of `SameCostsMarket{U}`, which contains the input data. Here
`U == eltype(t)` and we require `U <: Real`.

!!! note

    To improve the performance of the solvers, the schools are sorted in ascending order by `t`
    when the market is constructed. Therefore, `mkt = SameCostsMarket(f, t, h); mkt.t == t` is
    *not* necessarily true. The input data can be recovered using `mkt.t == t[mkt.perm]`; see
    [`SameCostsMarket(f, t, h)`](@ref).

To read in an instance of the varied-costs problem, use [`VariedCostsMarket(f, t, g, H)`](@ref):

```@example 1
mkt = VariedCostsMarket(
    [0.5, 0.1, 0.9, 0.7],   # admissions probabilities
    [12, 20, 1, 3],         # utility values
    [3, 4, 7, 2],           # application fees
    8                       # budget to spend on applications
)
```

The output is an instance of `VariedCostsMarket`. Both `SameCostsMarket` and `VariedCostsMarket` are subtypes
of the abstract type `Market`.

!!! note

    `VariedCostsMarket` supports only `t` with integer eltype. Because the objective function is linear in `t`,
    to work with float data, first multiply by the least common denominator.

If the application fees share a common divisor, the solvers will perform more effectively if you divide by it. For example,
fees of `[90, 80, 90, 70]` with budget `160` is equivalent to fees of `[9, 8, 9, 7]` with budget `16`.

OptimalApplication also provides the convenience functions [`Market(f, t, h)`](@ref) and
[`Market(f, t, g, H)`](@ref), which behave just like those above.

