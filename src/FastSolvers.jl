"""
    College(f, t, ft, omf)

Contains a college's admissions probability `f`, utility value `t`,
their product `ft`, and `1 - f = omf`. Used only by `applicationorder()`.
"""
struct College
    j::Int
    f::Float64
    t::Float64
    ft::Float64
    omf::Float64

    function College(j::Int, f::Float64, t::Real, ft::Float64, omf::Float64)
        return new(j, f, t, ft, omf)
    end

    function College(j::Int, f::Float64, t::Real, omf::Float64)
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
function applicationorder_list(mkt::SameCostsMarket)::Tuple{Vector{Int},Vector{Float64}}
    apporder = zeros(Int, mkt.h)
    v = zeros(mkt.h)

    mkt_list = College[College(j, mkt.f[j], mkt.t[j], mkt.ft[j], mkt.omf[j]) for j in 1:mkt.m]
    
    for j in 1:mkt.h
        c_k, k = findmax(mkt_list)
        v[j] = get(v, j - 1, 0) + c_k.ft
        apporder[j] = c_k.j

        for i in 1:k-1
            mkt_list[i] =
                College(
                    mkt_list[i].j,
                    mkt_list[i].f,
                    mkt_list[i].t * c_k.omf,
                    mkt_list[i].ft * c_k.omf,
                    mkt_list[i].omf
                )
        end
        for i in k+1:length(mkt_list)
            mkt_list[i] =
                College(
                    mkt_list[i].j,
                    mkt_list[i].f,
                    mkt_list[i].t - c_k.ft,
                    mkt_list[i].ft - mkt_list[i].f * c_k.ft,
                    mkt_list[i].omf
                )
        end
    
        deleteat!(mkt_list, k)
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
                # what we need. Equiv. to heap.drain!() in Rust
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
    X = Int64[]
    for j in mkt.m:-1:1
        if get(V, (j - 1, h), 0) < get(V, (j, h), 0)
            push!(X, j)
            h -= mkt.g[j]
        end
    end

    return X, V
end


const IntTypes = [Int8, Int16, Int32, Int64]



"""
    optimalportfolio_dynamicprogram(mkt::VariedCostsMarket)

Use the dynamic program on application costs to produce the optimal portfolio `X` and associated
value `v` for the market `mkt` with varying application costs. 
"""
function optimalportfolio_dynamicprogram(mkt::VariedCostsMarket)::Tuple{Vector{Int},Float64}
    V_dict = Dict{Tuple{Int, Int}, Float64}()

    function V(j, h)
        haskey(V_dict, (j, h)) && V_dict[(j, h)]

        if j == 0 || h == 0
            return 0.0
        elseif h < mkt.g[j]
            push!(V_dict, (j, h) => V(j - 1, h))
            return V_dict[(j, h)]
        else
            push!(V_dict, (j, h) => max(
                V(j - 1, h),
                mkt.omf[j] * V(j - 1, h - mkt.g[j]) + mkt.ft[j]
            ))
            return V_dict[(j, h)]
        end
    end

    h = mkt.H
    X = Int[]

    v = V(mkt.m, mkt.H)

    for j in mkt.m:-1:1
        if V(j - 1, h) < V(j, h)
            push!(X, j)
            h -= mkt.g[j]
        end
    end

    return X, v
end


"""
    optimalportfolio_fptas(mkt::VariedCostsMarket, ε)

Use the fully polynomial-time approximation scheme to produce a
`1-ε`-optimal portfolio for the market `mkt` with varying application costs. 
"""
function optimalportfolio_fptas(mkt::VariedCostsMarket, ε::Float64)::Tuple{Vector{Int},Float64}
    @assert 0 < ε < 1

    P = ceil(Int, log2(mkt.m^2 / (ε * sum(mkt.ft))))

    t = mkt.t .* 2^P
    ft = mkt.f .* t

    Ū = ceil(Int, sum(ft))

    infty = sum(mkt.g) + 1

    G_dict = Dict{Tuple{Int,Int},Int}()

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
                v_minus_Δ = floor(Int, clamp((v - ft[j]) / mkt.omf[j], -1, Ū))

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

    # Binary search
    v = 0
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

    X = Int[]

    for j in mkt.m:-1:1
        if G(j, v) < infty && G(j, v) < G(j - 1, v)
            push!(X, j)
            v = floor(Int, clamp((v - ft[j]) / mkt.omf[j], -1, Ū))
        end
    end

    return X, valuation(X, mkt)
end
