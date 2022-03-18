"""
    optimalportfolio_enumerate(mkt::SameCostsMarket)

Produce the optimal portfolio for the market `mkt` having identical application costs.
"""
function optimalportfolio_enumerate(mkt::SameCostsMarket)::Tuple{Vector{Int},Float64}
    X = zeros(Int64, mkt.h)
    v = 0.0

    for Y in multiset_combinations(1:mkt.m, mkt.h)
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
function optimalportfolio_enumerate(mkt::VariedCostsMarket)::Tuple{Vector{Int},Float64}
    let
        X = Set()
        v = 0.0
        for Y in combinations(1:mkt.m)
            if (w = valuation(Y, mkt)) > v && sum(mkt.g[Y]) ≤ mkt.H
                v = w
                X = copy(Y)
            end
        end

        return X, v
    end
end
