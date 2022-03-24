"""
    optimalportfolio_enumerate(mkt::SameCostsMarket)

Produce the optimal portfolio for the market `mkt` having identical application costs.
"""
function optimalportfolio_enumerate(mkt::SameCostsMarket{T})::Tuple{Vector{T},Float64} where T
    X = zeros(T, mkt.h)
    v = 0.0

    for Y in multiset_combinations(T(1):T(mkt.m), Int(mkt.h))
        if (w = valuation(Y, mkt)) > v
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
function optimalportfolio_enumerate(mkt::VariedCostsMarket{T})::Tuple{Vector{T},Float64} where T
    let
        X = Set{T}()
        v = 0.0
        for Y in combinations(T(1):mkt.m)
            if (w = valuation(Y, mkt)) > v && sum(mkt.g[Y]) ≤ mkt.H
                v = w
                X = copy(Y)
            end
        end

        return X, v
    end
end
