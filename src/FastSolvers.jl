"""
    College(f, t, ft, omf)

Contains a college's admissions probability `f`, utility value `t`,
their product `ft`, and `1 - f = omf`. Used only by `applicationorder()`.
"""
struct College
    j::Int16
    f::Float64
    t::Float64
    ft::Float64
    omf::Float64

    function College(j::Integer, f::Float64, t::Real, ft::Float64, omf::Float64)
        return new(Int16(j), f, t, ft, omf)
    end

    function College(j::Integer, f::Float64, t::Real, omf::Float64)
        return new(Int16(j), f, t, f * t, omf)
    end
end

# Overload isless() so that the heap is ordered by expected utility.
isless(c1::College, c2::College) = isless(c1.ft, c2.ft)

const dummy_college = College(0, 1.0, -1.0, -1.0, 0.0)


"""
    applicationorder_list(mkt::SameCostsMarket)

Produce the optimal order of application for the market `mkt` having identical
application costs and the corresponding portfolio valuations.
"""
function applicationorder_list(mkt::SameCostsMarket, verbose=false::Bool)::Tuple{Vector{Int},Vector{Float64}}
    apporder = zeros(Int, mkt.h)
    v = zeros(mkt.h)

    mkt_list = College[College(j, mkt.f[j], mkt.t[j], mkt.ft[j], mkt.omf[j]) for j in 1:mkt.m]
    
    c_best::College, idx_best::Int = findmax(mkt_list)
    for j in 1:mkt.h
        if verbose
            println("Iteration $j")
            println("  ft: ", [mkt_list[i].ft for i in 1:length(mkt_list)])
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
function applicationorder_heap(mkt::SameCostsMarket)::Tuple{Vector{Int},Vector{Float64}}
    apporder = zeros(Int, mkt.h)
    v = zeros(mkt.h)
    
    mkt_heap = BinaryMaxHeap{College}(collect(College(j, mkt.f[j], mkt.t[j], mkt.ft[j], mkt.omf[j]) for j in 1:mkt.m))
    
    for j in 1:mkt.h
        c_k = first(mkt_heap)
        v[j] = get(v, j - 1, 0) + c_k.ft
        apporder[j] = c_k.j
    
        mkt_heap = BinaryMaxHeap{College}(collect(
                College(
                    c.j,
                    c.f,
                    c.t < c_k.t ? c.t * c_k.omf : c.t - c_k.ft,
                    c.t < c_k.t ? c.ft * c_k.omf : c.ft - c.f * c_k.ft,
                    c.omf
                )
                # This heap.valtree field is not documented in DataStructures.jl,
                # but it contains the heap data in an arbitrary order, which is exactly 
                # what we need. Equiv. to heap.drain!() in Rust except doesn't mutate.
                for c in (d for d in mkt_heap.valtree if d.j != c_k.j)
            )
        )
    end

    return apporder, v
end


"""
    optimalportfolio_valuationtable(mkt::VariedCostsMarket)

Use the dynamic program on application costs to produce the optimal portfolio `X` and associated
valuation table `V` for the market `mkt` with varying application costs.
"""
function optimalportfolio_valuationtable(mkt::VariedCostsMarket)::Tuple{Vector{Int},Matrix{Float64}}
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
    X = Int[]
    for j in mkt.m:-1:1
        if get(V, (j - 1, h), 0) < get(V, (j, h), 0)
            push!(X, j)
            h -= mkt.g[j]
        end
    end

    return X, V
end


# Used by dynamic program below
function V_recursor!(V_dict::Dict{Tuple{Int16,Int},Float64}, j::Int16, h::Int, mkt::VariedCostsMarket)
    haskey(V_dict, (j, h)) && return V_dict[(j, h)]

    if j == 0 || h == 0
        return 0.0
    elseif h < mkt.g[j]
        push!(V_dict, (j, h) => V_recursor!(V_dict, j - Int16(1), h, mkt))
        return V_dict[(j, h)]
    else
        jmo = j - Int16(1)
        push!(V_dict, (j, h) => max(
            V_recursor!(V_dict, jmo, h, mkt),
            mkt.omf[j] * V_recursor!(V_dict, jmo, h - mkt.g[j], mkt) + mkt.ft[j]
        ))
        return V_dict[(j, h)]
    end
end



"""
    optimalportfolio_dynamicprogram(mkt::VariedCostsMarket)

Use the dynamic program on application costs to produce the optimal portfolio `X` and associated
value `v` for the market `mkt` with varying application costs. 
"""
function optimalportfolio_dynamicprogram(mkt::VariedCostsMarket)::Tuple{Vector{Int},Float64}
    V_dict = Dict{Tuple{Int16,Int},Float64}()
    v = V_recursor!(V_dict, Int16(mkt.m), mkt.H, mkt)
    
    h = mkt.H
    X = Int[]

    for j in Int16(mkt.m):Int16(-1):Int16(1)
        if V_recursor!(V_dict, j - Int16(1), h, mkt) < V_recursor!(V_dict, j, h, mkt)
            push!(X, j)
            h -= mkt.g[j]
        end
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
        P = ceil(Int, log2(mkt.m^2 / (ε * sum(mkt.ft))))
        t = mkt.t .* 2^P
        ft = mkt.f .* t
        Ū = ceil(Int, sum(ft))
        infty = sum(mkt.g) + 1

        return new(t, ft, Ū, infty)
    end
end


# Used by fptas below
function G_recursor!(
        G_dict::Dict{Tuple{Int16,Int},Int},
        j::Int16,
        v::Int, 
        mkt::VariedCostsMarket,
        sp::ScaleParams
    )::Int
    haskey(G_dict, (j, v)) && return G_dict[(j, v)]

    if v ≤ 0
        return 0
    elseif j == 0 || sp.t[j] < v || v ≥ sp.Ū
        return sp.infty
    else
        jmo = j - Int16(1)
        if mkt.f[j] < 1
            # Clamping prevents over/underflow: for any v≤0 or v>Ū the function
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


"""
    optimalportfolio_fptas(mkt::VariedCostsMarket, ε)

Use the fully polynomial-time approximation scheme to produce a
`1-ε`-optimal portfolio for the market `mkt` with varying application costs. 
"""
function optimalportfolio_fptas(mkt::VariedCostsMarket, ε::Float64)::Tuple{Vector{Int},Float64}
    sp = ScaleParams(mkt, ε)

    G_dict = Dict{Tuple{Int16,Int},Int}()

    # Binary search
    v = 0
    v_UB = sp.Ū

    m = Int16(mkt.m)

    while v + 1 < v_UB
        mid = (v + v_UB) ÷ 2

        if G_recursor!(G_dict, m, mid, mkt, sp) > mkt.H
            v_UB = mid
        else
            v = mid
            if G_recursor!(G_dict, m, v + 1, mkt, sp) > mkt.H
                @goto foundit
            end
        end
    end
    @label foundit

    X = Int[]

    for j in Int16(mkt.m):Int16(-1):Int16(1)
        if G_recursor!(G_dict, j, v, mkt, sp) < sp.infty && G_recursor!(G_dict, j, v, mkt, sp) < G_recursor!(G_dict, j - Int16(1), v, mkt, sp)
            push!(X, j)
            v = floor(Int, clamp((v - sp.ft[j]) / mkt.omf[j], -1, sp.Ū))
        end
    end

    return X, valuation(X, mkt)
end
