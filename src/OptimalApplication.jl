module OptimalApplication

using Combinatorics: combinations, multiset_combinations
using DataStructures
# For implementing the FPTAS for Ellis's problem
using FixedPointNumbers
import Base.isless


export
    valuation,
    applicationorder,
    optimalportfolio_dynamicprogram,
    optimalportfolio_enumerate,
    optimalportfolio_valuationtable,
    optimalportfolio_fptas


const IntTypes = [Int8, Int16, Int32, Int64]


"""
    College(f, t, ft)

Contains a college's admissions probability `f`, utility value `t`,
and their product `ft`.
"""
struct College
    f::Float64
    t::Real
    ft::Float64

    College(f::Float64, t::Real) = new(f, t, f * t)
end

# Overload isless() so that the heap is ordered by expected utility.
isless(c1::College, c2::College) = isless(c1.ft, c2.ft)


"""
    applicationorder(f, t, h=m; datastructure=:heap)

Produce the optimal order of application for the market defined by admissions probabilities `f`
and utility values `t`. Stores college data in a binary `:heap` or as a `:dict`. 
"""
function applicationorder(
    f::Vector{Float64},
    t::Vector{<:Real},
    h = nothing::Union{Int64,Nothing};
    datastructure = :heap::Symbol)::Tuple{Vector{Int64},Vector{Float64}}
    m = length(f)
    @assert m == length(t)
    @assert issorted(t)
    if isnothing(h)
        h = m
    else
        @assert 0 < h ≤ m
    end

    if datastructure == :heap
        mkt = MutableBinaryMaxHeap{College}()
        for j in 1:m
            push!(mkt, College(f[j], t[j]))
        end

        onheap = Set(1:m)
        apporder = zeros(Int, h)
        v = zeros(h)
        for j in 1:h
            c_k, k = top_with_handle(mkt)
            v[j] = get(v, j - 1, 0) + c_k.ft

            pop!(mkt)

            apporder[j] = k
            delete!(onheap, k)
            for i in onheap
                if i < k
                    update!(mkt, i, College(mkt[i].f, mkt[i].t * (1 - c_k.f)))
                else
                    update!(mkt, i, College(mkt[i].f, mkt[i].t - c_k.ft))
                end
            end
        end
    else
        mkt = Dict{Int64,College}()
        for j in 1:m
            push!(mkt, j => College(f[j], t[j]))
        end

        apporder = zeros(Int, h)
        v = zeros(h)
        for j in 1:h
            c_k, k = findmax(mkt)
            v[j] = get(v, j - 1, 0) + c_k.ft

            delete!(mkt, k)

            apporder[j] = k
            for i in keys(mkt)
                if i < k
                    mkt[i] = College(mkt[i].f, mkt[i].t * (1 - c_k.f))
                else
                    mkt[i] = College(mkt[i].f, mkt[i].t - c_k.ft)
                end
            end
        end
    end

    return apporder, v
end


"""
    valuation(X, f, t)

Returns the expected value of the portfolio `X`, a vector of school indices, on the
the market defined by admissions probabilities `f` and utility values `t`. 
"""
function valuation(
    X::Vector{Int64},
    f::Vector{<:Real},
    t::Vector{<:Real})::Float64
    @assert issorted(t)
    sort!(X)
    h = length(X)

    if h > 1
        res = 0.0
        cp = reverse(cumprod(1 .- reverse(f[X[2:end]])))

        for j in 1:h-1
            res += t[X[j]] * f[X[j]] * cp[j]
        end

        res += t[X[end]] * f[X[end]]

        return res
    else
        return f[X[1]] * t[X[1]]
    end
end


"""
    optimalportfolio_enumerate(f, t, h)

Produce the optimal portfolio of size `h` on the market defined by admissions probabilities `f`
and utility values `t`. Solves by enumeration. 
"""
function optimalportfolio_enumerate(
    f::Vector{<:Real},
    t::Vector{<:Real},
    h::Int64)::Tuple{Vector{Int64},Float64}
    m = length(t)
    X = zeros(Int64, h)
    v = 0.0

    for Y in multiset_combinations(1:m, h)
        if (w = valuation(Y, f, t)) > v
            v = w
            X[:] = Y
        end
    end

    return X, v
end


"""
    optimalportfolio_enumerate(f, t, g, H)

Produce the optimal portfolio of cost `H` on the market defined by admissions probabilities `f`,
utility values `t`, and application costs `g`. Solves by enumeration. 
"""
function optimalportfolio_enumerate(
    f::Vector{<:Real},
    t::Vector{<:Real},
    g::Vector{<:Real},
    H::Real)::Tuple{Vector{Int64},Float64}
    m = length(t)

    let
        X = Set()
        v = 0.0
        for Y in combinations(1:m)
            if (w = valuation(Y, f, t)) > v && sum(g[Y]) ≤ H
                v = w
                X = copy(Y)
            end
        end

        return X, v
    end
end


"""
    optimalportfolio_valuationtable(f, t, g, H)

Given admissions probabilities `f`, utility values `t`, application costs `g`, and
budget `H`, uses a dynamic program to produce the optimal portfolio `X` and associated
valuation table `V`.
"""
function optimalportfolio_valuationtable(
    f::Vector{Float64},
    t::Vector{<:Real},
    g::Vector{Int64},
    H::Int64)::Tuple{Vector{Int64},Matrix{Float64}}
    m = length(f)
    @assert m == length(t) == length(g)
    @assert issorted(t)
    @assert 0 < H ≤ sum(g)

    V = zeros(m, H)
    #     U = falses(m, H)
    for j in 1:m, h in 1:H
        if h < g[j]
            V[j, h] = get(V, (j - 1, h), 0)
        else
            V[j, h] = max(
                get(V, (j - 1, h), 0),
                (1 - f[j]) * get(V, (j - 1, h - g[j]), 0) + f[j] * t[j]
            )
        end
    end

    h = H
    X = Int64[]
    for j in m:-1:1
        if get(V, (j - 1, h), 0) < get(V, (j, h), 0)
            push!(X, j)
            h -= g[j]
        end
    end

    return X, V
end


"""
    optimalportfolio_dynamicprogram(f, t, g, H; memoize=true)

Given admissions probabilities `f`, utility values `t`, application costs `g`, and
budget `H`, uses a dynamic program to produce the optimal portfolio `X` and associated
value `v`. Set `memoize=false` to (unwisely) use blind recursion.
"""
function optimalportfolio_dynamicprogram(
    f::Vector{Float64},
    t::Vector{<:Real},
    g::Vector{Int64},
    H::Int64,
    memoize = true)::Tuple{Vector{Int64},Float64}
    m = length(f)
    @assert m == length(t) == length(g)
    @assert issorted(t)
    @assert 0 < H ≤ sum(g)

    T2 = IntTypes[findfirst(m < typemax(T) for T in IntTypes)]
    T3 = IntTypes[findfirst(H < typemax(T) for T in IntTypes)]

    if memoize
        V_dict = Dict{Tuple{T2,T3},Float64}()
    end

    function V(j, h)
        if memoize && haskey(V_dict, (j, h))
            return V_dict[(j, h)]
        end

        if j == 0 || h == 0
            return 0.0
        elseif h < g[j]
            if memoize
                push!(V_dict, (j, h) => V(j - 1, h))
                return V_dict[(j, h)]
            else
                return V(j - 1, h)
            end
        else
            if memoize
                push!(V_dict, (j, h) => max(
                    V(j - 1, h),
                    (1 - f[j]) * V(j - 1, h - g[j]) + f[j] * t[j]
                ))
                return V_dict[(j, h)]
            else
                return max(
                    V(j - 1, h),
                    (1 - f[j]) * V(j - 1, h - g[j]) + f[j] * t[j]
                )
            end
        end
    end

    h = H
    X = Int64[]

    if memoize
        v = V(m, H)
    end

    for j in m:-1:1
        if V(j - 1, h) < V(j, h)
            push!(X, j)
            h -= g[j]
        end
    end

    if memoize
        return X, v
    else
        return X, valuation(X, f, t)
    end
end


"""
    optimalportfolio_fptas(f, t, g, H)

Given admissions probabilities `f`, utility values `t`, application costs `g`, and
budget `H`, uses the fully polynomial-time approximation scheme to produce a
`1-ε`-optimal portfolio.
"""
function optimalportfolio_fptas(
    f::Vector{Float64},
    t::Vector{Int64},
    g::Vector{Int64},
    H::Int64,
    ε::Float64)::Tuple{Vector{Int64},Float64}
    m = length(f)
    sumg = sum(g)
    @assert m == length(t) == length(g)
    @assert issorted(t)
    @assert 0 < H ≤ sumg

    Ū = sum(f .* t)
    P = ceil(Int64, log2(m^2 / (ε * Ū)))
    infty = sumg + 1

    T1 = IntTypes[findfirst(Ū * 2^P < typemax(T) for T in IntTypes)]
    T2 = IntTypes[findfirst(m < typemax(T) for T in IntTypes)]
    T3 = IntTypes[findfirst(infty < typemax(T) for T in IntTypes)]

    FP = Fixed{T1,P}

    G_dict = Dict{Tuple{T2,FP},T3}()

    # Want to assert v::FP here but gives a strange error
    function G(j::Int, v::Real)
        haskey(G_dict, (j, v)) && return G_dict[(j, v)]

        if v ≤ 0
            return 0
        elseif j == 0 || t[j] < v
            push!(G_dict, (j, v) => infty)
            return infty # G_dict[(j, v)]
        else
            if f[j] < 1
                v_minus_Δ = FP(v - f[j] * (t[j] - v) / (1 - f[j]))
                push!(G_dict, (j, v) => min(
                    G(j - 1, v),
                    g[j] + G(j - 1, v_minus_Δ)
                ))
                return G_dict[(j, v)]
            else
                push!(G_dict, (j, v) => min(G(j - 1, v), g[j]))
                return G_dict[(j, v)]
            end
        end
    end

    # Binary search below is equivalent to the following linear search:
    # v = findlast(g -> g ≤ H, G(m, w) for w in eps(FP):eps(FP):FP(Ū)) * eps(FP)

    v = eps(FP)
    v_UB = FP(Ū)
    while v + eps(FP) < v_UB
        mid = (v + v_UB) / 2
        if G(m, mid) > H
            v_UB = mid
        else
            v = mid
            if G(m, v + eps(FP)) > H
                @goto foundit
            end
        end
    end
    @label foundit

    X = Int64[]

    for j in m:-1:1
        if G(j, v) < infty && G(j, v) < G(j - 1, v)
            push!(X, j)
            v = FP(v - f[j] * (t[j] - v) / (1 - f[j]))
        end
    end

    return X, valuation(X, f, t)
end


end