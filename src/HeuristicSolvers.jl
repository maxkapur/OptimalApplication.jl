function optimalportfolio_greedy_nopermute(mkt::VariedCostsMarket)::Vector{Int}
    priority_order = [j for j = 1:mkt.m if mkt.g[j] ≤ mkt.H]
    sort!(priority_order, by = function (j)
        mkt.f[j] * mkt.t[j] / mkt.g[j]
    end, rev = true)

    X = Int[]
    H = mkt.H
    for j in priority_order
        if mkt.g[j] ≤ H
            push!(X, j)
            H -= mkt.g[j]
        else
            break
        end
    end

    return X
end


"""
    optimalportfolio_greedy(mkt::VariedCostsMarket) -> X, v

Use a greedy heuristic that adds schools in decreasing order of `mkt.f .* mkt.t ./ mkt.g`
to compute a heuristically optimal portfolio for the `VariedCostsMarket`
defined by `mkt`.
"""
function optimalportfolio_greedy(mkt::VariedCostsMarket)::Tuple{Vector{Int},Float64}
    X = optimalportfolio_greedy_nopermute(mkt)
    sort!(X)
    return mkt.perm[X], valuation_nopermute_sorted(X, mkt)
end


"""
    neighbor(X, mkt) -> X_neighbor

Generate a random neighbor of `X` for the market `mkt`. If `X` is feasible,
preserve feasibility; if not, move to a feasible solution. `X` is assumed
to refer to the sorted indices of the market. 
"""
function neighbor(X::AbstractVector{<:Integer}, mkt::VariedCostsMarket)::Vector{Int}
    X_neighbor = copy(X)
    Y = setdiff(1:mkt.m, X)
    shuffle!(X_neighbor)
    shuffle!(Y)

    H = sum(mkt.g[X])

    while H ≤ mkt.H && !isempty(Y)
        j = pop!(Y)
        push!(X_neighbor, j)
        H += mkt.g[j]
    end

    while H > mkt.H
        j = popfirst!(X_neighbor)
        H -= mkt.g[j]
    end

    return X_neighbor
end

"""
    optimalportfolio_simulatedannealing(
        mkt::VariedCostsMarket;
        temp::Float64=0.25,
        red::Float64=0.0625,
        nit::Integer=500,
        verbose::Bool=false
    ) -> X, v

Use a simulated annealing procedure to compute a heuristically optimal portfolio
for the `VariedCostsMarket` defined by `mkt`. `temp` is the initial temperature,
`red` is the reduction factor, and `nit` is the total number of iterations.
"""
function optimalportfolio_simulatedannealing(
    mkt::VariedCostsMarket;
    temp::Float64 = 0.25,
    red::Float64 = 0.0625,
    nit::Integer = 500,
    verbose::Bool = false,
)::Tuple{Vector{Int},Float64}
    X = optimalportfolio_greedy_nopermute(mkt)
    sort!(X)
    v = valuation_nopermute_sorted(X, mkt)

    X_best, v_best = X, v

    for i = 1:nit
        X_neighbor = neighbor(X, mkt)
        sort!(X_neighbor)
        v_neighbor = valuation_nopermute_sorted(X_neighbor, mkt)

        if verbose
            println("Iteration $i, v_best = $v_best, v_neighbor = $v_neighbor")
        end

        if (delta = v_neighbor - v) ≥ 0
            X, v = X_neighbor, v_neighbor
            if v > v_best
                X_best, v_best = X, v
            end
        elseif rand() < exp(delta / temp)
            verbose && @show exp(delta / temp)
            X, v = X_neighbor, v_neighbor
        end

        temp *= red
    end

    return mkt.perm[X_best], v_best
end
