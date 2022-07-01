# Exact form of `G_recursor!` for use in `optimalportfolio_dynamicprogram_slow`.
# No memoization, since the frequency of index collisions is very low.
@inbounds function G_recursor_exact(
    j::Int,
    v::Float64,
    mkt::VariedCostsMarket,
    v_static_UB::Float64
)::Int
    if v ≤ 0
        return 0
    elseif j == 0 || mkt.t[j] < v || v ≥ v_static_UB
        return mkt.H + 1 # = infty, because infeasible
    else
        if mkt.f[j] < 1
            v_minus_Δ = (v - mkt.f[j] * mkt.t[j]) / (1 - mkt.f[j])

            return min(
                G_recursor_exact(j - 1, v, mkt, v_static_UB),
                mkt.g[j] + G_recursor_exact(j - 1, v_minus_Δ, mkt, v_static_UB)
            )
        else
            return min(G_recursor_exact(j - 1, v, mkt, v_static_UB), mkt.g[j])
        end
    end
end

"""
    optimalportfolio_dynamicprogram_slow(mkt, δ) -> X, v

Uses a dynamic program to produce the optimal portfolio `X`
with valuation `v`, for the [`VariedCostsMarket`](@ref) defined by `mkt`. 
The valuation of `X` differs from that of the global optimum by no more than `δ`.

```julia-repl
julia> mkt = VariedCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], [2, 4, 2, 5, 1], 8);

julia> optimalportfolio_dynamicprogram_slow(mkt, 1e-6)
([3, 5, 2], 3.24)
```

This dynamic program iterates on portfolio valuations and is an "exact" form
of the approximation algorithm implemented in this package as
[`optimalportfolio_fptas`](@ref). For the vast majority of problem instances,
[`optimalportfolio_dynamicprogram`](@ref) is a faster way to obtain an exact
solution, but this one technically can be used to obtain the global optimum by
setting `δ = eps()`.
"""
function optimalportfolio_dynamicprogram_slow(
    mkt::VariedCostsMarket,
    δ::Float64
)::Tuple{Vector{Int},Float64}
    v_static_UB = valuation(1:mkt.m, mkt)

    # Binary search for max{ w :  G_recursor_exact(mkt.m, w, mkt), sp ≤ H }
    v = 0.0
    v_current_UB = v_static_UB
    while v + δ < v_current_UB
        mid = (v + v_current_UB) / 2

        if G_recursor_exact(mkt.m, mid, mkt, v_static_UB) > mkt.H
            v_current_UB = mid
        else
            v = mid
        end
    end

    X = Int[]
    # Technically memoization would improve the speed of this loop,
    # but the binary search above takes so much longer that it's not worth it
    @inbounds for j in reverse(1:mkt.m)
        if G_recursor_exact(j, v, mkt, v_static_UB) < G_recursor_exact(j - 1, v, mkt, v_static_UB)
            push!(X, j)
            v = (v - mkt.f[j] * mkt.t[j]) / (1 - mkt.f[j])
        end
    end

    return mkt.perm[X], valuation_nopermute_sorted(reverse(X), mkt)
end
