module OptimalApplication

using Combinatorics: combinations, multiset_combinations
using DataStructures
import Base.isless
import Base.hash


export
    Market,
    SameCostsMarket,
    VariedCostsMarket,
    valuation,
    applicationorder,
    optimalportfolio_dynamicprogram,
    optimalportfolio_enumerate,
    optimalportfolio_valuationtable,
    optimalportfolio_fptas,
    optimalportfolio_branchbound


const IntTypes = [Int8, Int16, Int32, Int64]


"""
    College(f, t, ft)

Contains a college's admissions probability `f`, utility value `t`,
and their product `ft`.
"""
struct College
    f::Real
    t::Real
    ft::Real
    omf::Real

    function College(f::Real, t::Real, ft::Real, omf::Real)
        return new(f, t, ft, omf)
    end

    function College(f::Real, t::Real, omf::Real)
        return new(f, t, f * t, omf)
    end
end

# Overload isless() so that the heap is ordered by expected utility.
isless(c1::College, c2::College) = isless(c1.ft, c2.ft)


function iscoherentmarket(f::Vector{<:Real}, t::Vector{<:Real})
    @assert length(f) == length(t)
    @assert all(0 .< f .≤ 1)
    @assert issorted(t)
end


function iscoherentmarket(f::Vector{<:Real}, t::Vector{<:Real}, g::Vector{<:Real})
    @assert length(f) == length(t)
    @assert length(f) == length(g)
    @assert all(0 .< f .≤ 1)
    @assert issorted(t)
end


function isnontrivialmarket(f::Vector{<:Real}, t::Vector{<:Real}, h::Int)
    iscoherentmarket(f, t)
    @assert 0 < h ≤ length(t)
end


function isnontrivialmarket(
    f::Vector{<:Real},
    t::Vector{<:Real},
    g::Vector{<:Real},
    H::Real)
    iscoherentmarket(f, t, g)
    @assert 0 < H ≤ sum(g)
    @assert all(0 .< g .≤ H)
end


"""
Contains information about a college application market with identical
application costs.
"""
struct SameCostsMarket
    m::Integer
    f::Vector{<:Real}
    t::Vector{<:Real}
    h::Integer
    ft::Vector{<:Real}      # = f .* t
    omf::Vector{<:Real}     # = 1 .- f
    perm::Vector{<:Integer}

    function SameCostsMarket(
        f::Vector{<:Real},
        t::Vector{<:Real},
        h::Integer)
        isnontrivialmarket(f, t, h)
        m = length(f)

        # We are current asserting that t is sorted; a placeholder for later work
        perm = sortperm(t)

        return new(m, f, t, h, f .* t, 1 .- f, perm)
    end
end


"""
Contains information about a college application market with varying
application costs.
"""
struct VariedCostsMarket
    m::Integer
    f::Vector{<:Real}
    t::Vector{<:Real}
    g::Vector{<:Real}
    H::Real
    ft::Vector{<:Real}      # = f .* t
    omf::Vector{<:Real}     # = 1 .- f
    perm::Vector{<:Integer}

    function VariedCostsMarket(
        f::Vector{<:Real},
        t::Vector{<:Real},
        g::Vector{<:Real},
        H::Real)
        isnontrivialmarket(f, t, g, H)
        m = length(f)
        perm = sortperm(t)

        return new(m, f, t, g, H, f .* t, 1 .- f, perm)
    end
end


function Market(f::Vector{<:Real}, t::Vector{<:Real}, h::Integer)
    return SameCostsMarket(f, t, h)
end


function Market(f::Vector{<:Real}, t::Vector{<:Real}, g::Vector{<:Real}, H::Real)
    return VariedCostsMarket(f, t, g, H)
end


function valuation(X::Vector{<:Integer}, mkt::Union{SameCostsMarket,VariedCostsMarket})
    isempty(X) && return 0.0

    sort!(X)
    h = length(X)

    if h > 1
        res = 0.0
        cp = reverse(cumprod(reverse(mkt.omf[X[2:end]])))

        for j in 1:h-1
            res += mkt.ft[X[j]] * cp[j]
        end

        res += mkt.ft[X[end]]

        return res
    else
        return mkt.ft[X[1]]
    end
end


"""
    applicationorder(f, t, h=m; datastructure=:heap)

Produce the optimal order of application for the market defined by admissions probabilities `f`
and utility values `t`. Stores college data in a binary `:heap` or as a `:dict`. 
"""
function applicationorder(
    mkt::SameCostsMarket;
    datastructure = :heap::Symbol)::Tuple{Vector{Int64},Vector{Float64}}

    apporder = zeros(Int, mkt.h)
    v = zeros(mkt.h)
    if datastructure == :heap
        mkt_heap = MutableBinaryMaxHeap{College}()
        for j in 1:mkt.m
            push!(mkt_heap, College(mkt.f[j], mkt.t[j], mkt.ft[j], mkt.omf[j]))
        end

        onheap = Set(1:mkt.m)
        for j in 1:mkt.h
            c_k, k = top_with_handle(mkt_heap)
            v[j] = get(v, j - 1, 0) + c_k.ft

            pop!(mkt_heap)

            apporder[j] = k
            delete!(onheap, k)
            for i in onheap
                if i < k
                    update!(mkt_heap, i, College(mkt_heap[i].f, mkt_heap[i].t * c_k.omf, mkt_heap[i].omf))
                else
                    update!(mkt_heap, i, College(mkt_heap[i].f, mkt_heap[i].t - c_k.ft, mkt_heap[i].omf))
                end
            end
        end
    else
        mkt_dict = Dict{Int64,College}()
        for j in 1:mkt.m
            push!(mkt_dict, j => College(mkt.f[j], mkt.t[j], mkt.ft[j], mkt.omf[j]))
        end

        for j in 1:mkt.h
            c_k, k = findmax(mkt_dict)
            v[j] = get(v, j - 1, 0) + c_k.ft

            delete!(mkt_dict, k)

            apporder[j] = k
            for i in keys(mkt_dict)
                if i < k
                    mkt_dict[i] = College(mkt_dict[i].f, mkt_dict[i].t * c_k.omf, mkt_dict[i].omf)
                else
                    mkt_dict[i] = College(mkt_dict[i].f, mkt_dict[i].t - c_k.ft, mkt_dict[i].omf)
                end
            end
        end
    end

    return apporder, v
end


"""
    optimalportfolio_enumerate(f, t, h)

Produce the optimal portfolio of size `h` on the market defined by admissions probabilities `f`
and utility values `t`. Solves by enumeration. 
"""
function optimalportfolio_enumerate(mkt::SameCostsMarket)::Tuple{Vector{Int64},Float64}
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
    optimalportfolio_enumerate(f, t, g, H)

Produce the optimal portfolio of cost `H` on the market defined by admissions probabilities `f`,
utility values `t`, and application costs `g`. Solves by enumeration. 
"""
function optimalportfolio_enumerate(mkt::VariedCostsMarket)::Tuple{Vector{Int64},Float64}
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


# function optimalportfolio_enumerate(mkt::SameCostsMarket)
#     return optimalportfolio_enumerate(mkt.f, mkt.t, mkt.h)
# end


# function optimalportfolio_enumerate(mkt::VariedCostsMarket)
#     return optimalportfolio_enumerate(mkt.f, mkt.t, mkt.g, mkt.H)
# end


"""
    optimalportfolio_valuationtable(f, t, g, H)

Given admissions probabilities `f`, utility values `t`, application costs `g`, and
budget `H`, uses a dynamic program to produce the optimal portfolio `X` and associated
valuation table `V`.
"""
function optimalportfolio_valuationtable(mkt::VariedCostsMarket)::Tuple{Vector{Int64},Matrix{Float64}}
    @assert eltype(mkt.g) <: Integer
    @assert typeof(mkt.H) <: Integer

    V = zeros(mkt.m, mkt.H)
    for j in 1:mkt.m, h in 1:mkt.H
        if h < mkt.g[j]
            V[j, h] = get(V, (j - 1, h), 0)
        else
            V[j, h] = max(
                get(V, (j - 1, h), 0),
                mkt.omf[j] * get(V, (j - 1, h - mkt.g[j]), 0) + mkt.ft[j]
            )
        end
    end

    h = mkt.H
    X = Int64[]
    for j in mkt.m:-1:1
        if get(V, (j - 1, h), 0) < get(V, (j, h), 0)
            push!(X, j)
            h -= mkt.g[j]
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
function optimalportfolio_dynamicprogram(mkt::VariedCostsMarket, memoize = true)::Tuple{Vector{Int64},Float64}
    @assert eltype(mkt.g) <: Integer
    @assert typeof(mkt.H) <: Integer

    T2 = IntTypes[findfirst(mkt.m < typemax(T) for T in IntTypes)]
    T3 = IntTypes[findfirst(mkt.H < typemax(T) for T in IntTypes)]

    if memoize
        V_dict = Dict{Tuple{T2,T3},Float64}()
    end

    function V(j, h)
        if memoize && haskey(V_dict, (j, h))
            return V_dict[(j, h)]
        end

        if j == 0 || h == 0
            return 0.0
        elseif h < mkt.g[j]
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
                    mkt.omf[j] * V(j - 1, h - mkt.g[j]) + mkt.ft[j]
                ))
                return V_dict[(j, h)]
            else
                return max(
                    V(j - 1, h),
                    mkt.omf[j] * V(j - 1, h - mkt.g[j]) + mkt.ft[j]
                )
            end
        end
    end

    h = mkt.H
    X = Int64[]

    if memoize
        v = V(mkt.m, mkt.H)
    end

    for j in mkt.m:-1:1
        if V(j - 1, h) < V(j, h)
            push!(X, j)
            h -= mkt.g[j]
        end
    end

    if memoize
        return X, v
    else
        return X, valuation(X, mkt)
    end
end


"""
    optimalportfolio_fptas(f, t, g, H)

Given admissions probabilities `f`, utility values `t`, application costs `g`, and
budget `H`, uses the fully polynomial-time approximation scheme to produce a
`1-ε`-optimal portfolio.
"""
function optimalportfolio_fptas(
    mkt::VariedCostsMarket,
    ε::Float64)::Tuple{Vector{Int64},Float64}
    @assert eltype(mkt.t) <: Integer
    @assert 0 < ε < 1

    P = ceil(Int, log2(mkt.m^2 / (ε * sum(mkt.ft))))

    t = mkt.t .* 2^P
    ft = mkt.f .* t

    Ū = ceil(Int, sum(ft))

    infty = sum(mkt.g) + 1

    vType = IntTypes[findfirst(Ū < typemax(T) for T in IntTypes)]
    jType = IntTypes[findfirst(mkt.m < typemax(T) for T in IntTypes)]
    gType = IntTypes[findfirst(infty < typemax(T) for T in IntTypes)]

    Ū = vType(Ū)
    G_dict = Dict{Tuple{jType,vType},gType}()

    function G(j::Integer, v::Integer)
        haskey(G_dict, (j, v)) && return G_dict[(j, v)]

        if v ≤ 0
            return 0
        elseif j == 0 || t[j] < v || v ≥ Ū
            return infty
        else
            if mkt.f[j] < 1
                # Clamping prevents over/underflow: for any v≤0 or v>Ū the function
                # is trivially defined, so recording any more extreme number is meaningless
                v_minus_Δ = floor(vType, clamp((v - ft[j]) / mkt.omf[j], -1, Ū))

                push!(G_dict, (j, v) => min(
                    G(j - 1, v),
                    mkt.g[j] + G(j - 1, v_minus_Δ)
                ))
                return G_dict[(j, v)]
            else
                push!(G_dict, (j, v) => min(G(j - 1, v), mkt.g[j]))
                return G_dict[(j, v)]
            end
        end
    end

    # Binary search below is equivalent to the following linear search:
    # v = findlast(g -> g ≤ H, G(m, w) for w in eps(FP):eps(FP):FP(Ū)) * eps(FP)

    v = vType(0)
    v_UB = Ū

    while v + 1 < v_UB
        mid = (v + v_UB) ÷ 2

        if G(mkt.m, mid) > mkt.H
            v_UB = mid
        else
            v = mid
            if G(mkt.m, v + 1) > mkt.H
                @goto foundit
            end
        end
    end
    @label foundit

    X = Int64[]

    for j in mkt.m:-1:1
        if G(j, v) < infty && G(j, v) < G(j - 1, v)
            push!(X, j)
            v = floor(vType, clamp((v - ft[j]) / mkt.omf[j], -1, Ū))
        end
    end

    return X, valuation(X, mkt)
end


include("VariedCostsBnB.jl")


end
