"""
    optimalportfolio_greedy(mkt::VariedCostsMarket) -> X, v

Use a greedy heuristic that adds schools in decreasing order of `mkt.ft ./ mkt.g`
to compute an approximately optimal portfolio for the market `mkt` with varying application costs.
"""
function optimalportfolio_greedy(mkt::VariedCostsMarket{T})::Tuple{Vector{Int},Float64} where {T}
    priority_order = filter(
        sortperm(mkt.ft ./ mkt.g, rev=true)
    ) do j
        mkt.g[j] ≤ mkt.H
    end

    X = T[]
    H = mkt.H
    for j in priority_order
        if mkt.g[j] ≤ H
            push!(X, j)
            H -= mkt.g[j]
        end

        if H == 0
            @goto full
        end
    end
    @label full

    # In the valuation we just use the identity permutation as invp
    # To prevent wasteful permuting and then invpermuting
    return mkt.perm[X], valuation(X, mkt; invp=oneunit(T):mkt.m)
end


"""
    neighbor(X, mkt) -> X_neighbor

Generate a random neighbor of `X` for the market `mkt`. If `X` is feasible,
preserve feasibility; if not, move to a feasible solution. `X` is assumed
to refer to the sorted indices of the market. 
"""
function neighbor(X::AbstractVector{T}, mkt::VariedCostsMarket{T})::Vector{T} where {T}
    X_neighbor = copy(X)
    Y = setdiff(oneunit(T):mkt.m, X)
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
        X0::Union{Nothing,Vector{T}}=nothing,
        temp::Union{Nothing,Float64}=nothing,
        nit::Integer=500,
        red::Float64=0.95,
        verbose::Bool=false
    ) -> X, v

Use a simulated annealing procedure to compute a heuristically optimal portfolio
for `mkt` and its valuation. 
"""
function optimalportfolio_simulatedannealing(
    mkt::VariedCostsMarket{T};
    X0::Union{Nothing,Vector{T}}=nothing,
    temp::Float64=0.25,
    nit::Integer=500,
    red::Float64=0.0625,
    verbose::Bool=false
)::Tuple{Vector{Int},Float64} where {T}
    if X0 === nothing
        # Default to the greedy solution. We code it separately so we can get
        # X in permuted order rather than original.
        X, v = begin
            priority_order = filter(
                sortperm(mkt.ft ./ mkt.g, rev=true)
            ) do j
                mkt.g[j] ≤ mkt.H
            end
    
            X = T[]
            H = mkt.H
            for j in priority_order
                if mkt.g[j] ≤ H
                    push!(X, j)
                    H -= mkt.g[j]
                else
                    @goto full
                end
            end
            @label full
    
            X, valuation(X, mkt; invp=oneunit(T):mkt.m)
        end
    else
        X, v = copy(X0), valuation(X0, mkt)
    end

    X_best, v_best = X, v

    for i in 1:nit
        X_neighbor = neighbor(X, mkt)
        v_neighbor = valuation(X_neighbor, mkt; invp=oneunit(T):mkt.m)

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
