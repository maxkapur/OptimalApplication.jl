"""
    College(f, t, ft, omf)

Contains a college's admissions probability `f`, utility value `t`,
their product `ft`, and `1 - f = omf`. Used only by `applicationorder()`.
"""
struct College{T<:Unsigned}
    j::T
    f::Float64
    t::Float64
    ft::Float64
    omf::Float64

    function College(j::T, f::Float64, t::Real, ft::Float64, omf::Float64) where T<:Unsigned
        return new{T}(j, f, t, ft, omf)
    end

    function College{T}(j::T, f::Float64, t::Real, omf::Float64) where T<:Unsigned
        return new(j, f, t, f * t, omf)
    end
end

# Overload isless() so that the heap is ordered by expected utility.
isless(c1::College, c2::College) = isless(c1.ft, c2.ft)


"""
    applicationorder_list(mkt::SameCostsMarket)

Produce the optimal order of application for the market `mkt` having identical
application costs and the corresponding portfolio valuations.
"""
function applicationorder_list(mkt::SameCostsMarket{T}; verbose=false::Bool)::Tuple{Vector{Int},Vector{Float64}} where T
    apporder = zeros(T, mkt.h)
    v = zeros(mkt.h)

    mkt_list = College{T}[College(T(j), mkt.f[j], mkt.t[j], mkt.ft[j], mkt.omf[j]) for j in T(1):mkt.m]

    dummy_college = College(T(0), 1.0, -1.0, -1.0, 0.0)

    c_best::College{T}, idx_best::Int = findmax(mkt_list)
    @inbounds for j in T(1):mkt.h
        if verbose
            println("Iteration $j")
            println("  ft: ", [mkt_list[i].ft for i in T(1):length(mkt_list)])
            println("  Add school $(c_best.j)")
        end

        v[j] = get(v, j - 1, 0) + c_best.ft
        apporder[j] = c_best.j
    
        next_c_best::College = dummy_college
        next_idx_best::Int = 0
    
        for i in 1:idx_best-1
            mkt_list[i] =
                College(
                    mkt_list[i].j,
                    mkt_list[i].f,
                    mkt_list[i].t * c_best.omf,
                    mkt_list[i].ft * c_best.omf,
                    mkt_list[i].omf
                )
    
            if isless(next_c_best, mkt_list[i])
                next_c_best = mkt_list[i]
                next_idx_best = i
            end
        end
    
        for i in idx_best+1:length(mkt_list)
            mkt_list[i-1] =
                College(
                    mkt_list[i].j,
                    mkt_list[i].f,
                    mkt_list[i].t - c_best.ft,
                    mkt_list[i].ft - mkt_list[i].f * c_best.ft,
                    mkt_list[i].omf
                )
    
            if isless(next_c_best, mkt_list[i-1])
                next_c_best = mkt_list[i-1]
                next_idx_best = i - 1
            end
        end
    
        c_best, idx_best = next_c_best, next_idx_best
    
        pop!(mkt_list)
    end


    return apporder, v
end


"""
    applicationorder_heap(mkt::SameCostsMarket)

Produce the optimal order of application for the market `mkt` having identical
application costs and the corresponding portfolio valuations.
"""
function applicationorder_heap(mkt::SameCostsMarket{T})::Tuple{Vector{Int},Vector{Float64}} where {T}
    apporder = zeros(T, mkt.h)
    v = zeros(mkt.h)

    mkt_heap = BinaryMaxHeap{College{T}}(collect(College(T(j), mkt.f[j], mkt.t[j], mkt.ft[j], mkt.omf[j]) for j in T(1):mkt.m))

    @inbounds for j in T(1):mkt.h
        c_k = first(mkt_heap)
        v[j] = get(v, j - 1, 0) + c_k.ft
        apporder[j] = c_k.j
    
        mkt_heap = BinaryMaxHeap{College{T}}(
            map(filter(c -> c.j != c_k.j, mkt_heap.valtree)) do c
                College(
                    c.j,
                    c.f,
                    c.t < c_k.t ? c.t * c_k.omf : c.t - c_k.ft,
                    c.t < c_k.t ? c.ft * c_k.omf : c.ft - c.f * c_k.ft,
                    c.omf
                )
            end
        )
    end

    return apporder, v
end


"""
    optimalportfolio_valuationtable(mkt::VariedCostsMarket)

Use the dynamic program on application costs to produce the optimal portfolio `X` and associated
valuation table `V` for the market `mkt` with varying application costs.
"""
function optimalportfolio_valuationtable(mkt::VariedCostsMarket{T})::Tuple{Vector{Int},Matrix{Float64}} where {T}
    V = zeros(mkt.m, mkt.H)
    @inbounds for j in T(1):mkt.m, h in 1:mkt.H
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
    X = T[]
    for j in reverse(1:mkt.m)
        if get(V, (j - 1, h), 0) < get(V, (j, h), 0)
            push!(X, j)
            h -= mkt.g[j]
        end
    end

    return X, V
end


# Used by dynamic program below
@inbounds function V_recursor!(V_dict::Dict{Tuple{T,Int},Float64}, j::T, h::Int, mkt::VariedCostsMarket{T})::Float64 where T
    return get(V_dict, (j, h)) do
        if j == 0 || h == 0
            return 0.0
        elseif h < mkt.g[j]
            push!(V_dict, (j, h) => V_recursor!(V_dict, j - T(1), h, mkt))
            return V_dict[(j, h)]
        else
            jmo = j - T(1)
            push!(V_dict, (j, h) => max(
                V_recursor!(V_dict, jmo, h, mkt),
                mkt.omf[j] * V_recursor!(V_dict, jmo, h - mkt.g[j], mkt) + mkt.ft[j]
            ))
            return V_dict[(j, h)]
        end
    end
end



"""
    optimalportfolio_dynamicprogram(mkt::VariedCostsMarket)

Use the dynamic program on application costs to produce the optimal portfolio `X` and associated
value `v` for the market `mkt` with varying application costs. 
"""
function optimalportfolio_dynamicprogram(mkt::VariedCostsMarket{T}; verbose=false::Bool)::Tuple{Vector{Int},Float64} where {T}
    V_dict = Dict{Tuple{T,Int},Float64}()
    sizehint!(V_dict, mkt.m * mkt.m ÷ 2)

    v = V_recursor!(V_dict, mkt.m, mkt.H, mkt)

    h = mkt.H
    X = T[]

    for j in reverse(1:mkt.m)
        if V_recursor!(V_dict, T(j - 1), h, mkt) < V_recursor!(V_dict, T(j), h, mkt)
            push!(X, j)
            h -= mkt.g[j]
        end
    end

    if verbose
        V_table = Array{Union{Missing,Float64}}(missing, mkt.m, mkt.H)
        for (j, h) in keys(V_dict)
            if 0 < j ≤ mkt.m && 0 < h ≤ mkt.H
                V_table[j, h] = V_dict[(j, h)]
            end
        end

        display(V_table)
    end

    return X, v
end


# Used by fptas below
struct ScaleParams
    t::Vector{Int}
    ft::Vector{Float64}
    Ū::Int
    infty::Int

    function ScaleParams(mkt::VariedCostsMarket, ε::Float64)
        @assert 0 < ε < 1
        P = ceil(Int, log2(Int(mkt.m)^2 / (ε * sum(mkt.ft))))
        t = mkt.t .* 2^P
        ft = mkt.f .* t
        Ū = ceil(Int, sum(ft))
        infty = sum(mkt.g) + 1

        return new(t, ft, Ū, infty)
    end
end


# Used by fptas below
@inbounds function G_recursor!(
        G_dict::Dict{Tuple{T,Int},Int},
        j::T,
        v::Int, 
        mkt::VariedCostsMarket{T},
        sp::ScaleParams
    )::Int where T<:Unsigned
    return get(G_dict, (j, v)) do
        if v ≤ 0
            return 0
        elseif j == 0 || sp.t[j] < v || v ≥ sp.Ū
            return sp.infty
        else
            jmo = j - T(1)
            if mkt.f[j] < 1
                # Clamping prevents over/underflow: for any v<0 or v≥Ū the function
                # is trivially defined, so recording any more extreme number is meaningless
                v_minus_Δ = floor(Int, clamp((v - sp.ft[j]) / mkt.omf[j], -1, sp.Ū))

                push!(G_dict, (j, v) => min(
                    G_recursor!(G_dict, jmo, v, mkt, sp),
                    mkt.g[j] + G_recursor!(G_dict, jmo, v_minus_Δ, mkt, sp)
                ))
                return G_dict[(j, v)]
            else
                push!(G_dict, (j, v) => min(G_recursor!(G_dict, jmo, v, mkt, sp), mkt.g[j]))
                return G_dict[(j, v)]
            end
        end
    end
end


"""
    optimalportfolio_fptas(mkt::VariedCostsMarket, ε)

Use the fully polynomial-time approximation scheme to produce a
`1-ε`-optimal portfolio for the market `mkt` with varying application costs. 
"""
function optimalportfolio_fptas(mkt::VariedCostsMarket{T}, ε::Float64; verbose=false::Bool)::Tuple{Vector{Int},Float64} where T<:Unsigned
    sp = ScaleParams(mkt, ε)

    G_dict = Dict{Tuple{T,Int},Int}()
    sizehint!(G_dict, mkt.m * mkt.m ÷ 2)

    # Binary search
    v = 0
    v_UB = sp.Ū

    while v + 1 < v_UB
        mid = (v + v_UB) ÷ 2

        if G_recursor!(G_dict, mkt.m, mid, mkt, sp) > mkt.H
            v_UB = mid
        else
            v = mid
            if G_recursor!(G_dict, mkt.m, v + 1, mkt, sp) > mkt.H
                @goto foundit
            end
        end
    end
    @label foundit

    X = T[]

    for j in reverse(1:mkt.m)
        # G_recursor!(G_dict, T(j), v, mkt, sp) < sp.infty &&
        if G_recursor!(G_dict, T(j), v, mkt, sp) < G_recursor!(G_dict, T(j - 1), v, mkt, sp)
            push!(X, j)
            v = floor(Int, clamp((v - sp.ft[j]) / mkt.omf[j], -1, sp.Ū))
        end
    end

    if verbose
        G_table = Array{Union{Missing,Int}}(missing, mkt.m, v_UB)
        for (j, v) in keys(G_dict)
            if 0 < j ≤ mkt.m && 0 < v ≤ v_UB
                G_table[j, v] = G_dict[(j, v)]
            end
        end
    
        display(G_table)
    end


    return X, valuation(X, mkt)
end
