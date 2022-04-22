"""
    optimalportfolio_enumerate(mkt::SameCostsMarket) -> X, v

Produce the optimal portfolio for [`SameCostsMarket`](@ref) defined by `mkt`. Very slow;
you probably should use [`applicationorder_list`](@ref) or [`applicationorder_heap`](@ref) instead.

```julia-repl
julia> mkt = SameCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], 2);

julia> optimalportfolio_enumerate(mkt)
([2, 3], 2.7)
```
"""
function optimalportfolio_enumerate(mkt::SameCostsMarket)::Tuple{Vector{Int},Float64}
    mkt.m ≥ 21 && @warn "Enumeration is slow for large markets"

    X = zeros(Int, mkt.h)
    v = 0.0

    for Y in multiset_combinations(1:mkt.m, mkt.h)
        if (w = valuation_nopermute(Y, mkt)) > v
            v = w
            X[:] = Y
        end
    end

    return mkt.perm[X], v
end


"""
    optimalportfolio_enumerate(mkt::VariedCostsMarket) -> X, v

Produce the optimal portfolio for [`VariedCostsMarket`](@ref) defined by `mkt`. Very slow;
you probably should use [`optimalportfolio_dynamicprogram`](@ref) or [`optimalportfolio_branchbound`](@ref)
or [`optimalportfolio_simulatedannealing`](@ref) instead.

```julia-repl
julia> mkt = VariedCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], [2, 4, 2, 5, 1], 8);

julia> optimalportfolio_enumerate(mkt)
([2, 5, 3], 3.24)
```
"""
function optimalportfolio_enumerate(mkt::VariedCostsMarket)::Tuple{Vector{Int},Float64}
    mkt.m ≥ 21 && @warn "Enumeration is slow for large markets"

    X = Int[]
    v = 0.0

    for Y in combinations(1:mkt.m)
        if (w = valuation_nopermute(Y, mkt)) > v && sum(mkt.g[Y]) ≤ mkt.H
            v = w
            X = Y
        end
    end

    return mkt.perm[X], v
end
