"""
    optimalportfolio_enumerate(mkt::SameCostsMarket)

Produce the optimal portfolio for the market `mkt` having identical application costs.
"""
function optimalportfolio_enumerate(mkt::SameCostsMarket{T,U})::Tuple{Vector{Int},Float64} where {T,U}
    mkt.m ≥ 21 && @warn "Enumeration is slow for large markets"

    X = zeros(T, mkt.h)
    v = 0.0

    invp = invperm(mkt.perm)
    for Y in multiset_combinations(T(1):T(mkt.m), Int(mkt.h))
        if (w = valuation(Y, mkt; invp=invp)) > v
            v = w
            X[:] = Y
        end
    end

    return X, v
end


"""
    optimalportfolio_enumerate(mkt::VariedCostsMarket)

Produce the optimal portfolio for the market `mkt` having varying application costs.
"""
function optimalportfolio_enumerate(mkt::VariedCostsMarket{T})::Tuple{Vector{Int},Float64} where {T}
    mkt.m ≥ 21 && @warn "Enumeration is slow for large markets"

    let
        X = Set{T}()
        v = 0.0

        invp = invperm(mkt.perm)
        for Y in combinations(T(1):mkt.m)
            if (w = valuation(Y, mkt, invp=invp)) > v && sum(mkt.g[invp[Y]]) ≤ mkt.H
                v = w
                X = copy(Y)
            end
        end

        return X, v
    end
end
