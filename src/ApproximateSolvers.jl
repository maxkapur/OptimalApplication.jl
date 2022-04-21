# Used by fptas below
struct ScaleParams
    t::Vector{Int}
    ft::Vector{Float64}
    Ū::Int
    infty::Int

    function ScaleParams(mkt::VariedCostsMarket, ε::Float64)
        @assert 0 < ε < 1
        P = max(0,
            ceil(Int, log2(Int(mkt.m)^2 / (ε * sum(mkt.ft))))
        )
        t = mkt.t .* 2^P
        ft = mkt.f .* t
        Ū = ceil(Int, sum(ft))
        infty = sum(mkt.g) + 1

        return new(t, ft, Ū, infty)
    end
end


# Used by fptas below
@inbounds function G_recursor!(
    G_dict::Dict{Tuple{Int,Int},Int},
    j::Int,
    v::Int,
    mkt::VariedCostsMarket,
    sp::ScaleParams
)::Int
    return get(G_dict, (j, v)) do
        if v ≤ 0
            return 0
        elseif j == 0 || sp.t[j] < v || v ≥ sp.Ū
            return sp.infty
        else
            jmo = j - 1
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


# Function I stole from Base.Sort that does midpoint using bitshift.
# Performance-optimized but safe only if `lo <= hi`.
midpoint(lo::T, hi::T) where {T<:Integer} = lo + ((hi - lo) >>> 0x01)


"""
    optimalportfolio_fptas(mkt, ε; verbose=false) -> X, v

Use the fully polynomial-time approximation scheme to produce a `1-ε`-optimal portfolio `X`,
with valuation `v`, for the [`VariedCostsMarket`](@ref) defined by `mkt`.
"""
function optimalportfolio_fptas(
    mkt::VariedCostsMarket,
    ε::Float64;
verbose::Bool=false
)::Tuple{Vector{Int},Float64}
    sp = ScaleParams(mkt, ε)

    G_dict = Dict{Tuple{Int,Int},Int}()
    sizehint!(G_dict, mkt.m * mkt.m ÷ 2)

    # Binary search for max{ w :  G_recursor!(G_dict, mkt.m, w, mkt, sp ≤ H }

    # Should be able to use something like this with a proper by kw in searchsortedlast
    # vs = 0:sp.Ū
    # searchsortedlast(vs)


    v = 0
    v_UB = sp.Ū

    @inbounds while v + 1 < v_UB
        # mid = (v + v_UB) ÷ 2
        mid = midpoint(v, v_UB)

        if G_recursor!(G_dict, mkt.m, mid, mkt, sp) > mkt.H
            v_UB = mid
        else
            v = mid
        end
    end

    X = Int[]

    for j in reverse(1:mkt.m)
        # G_recursor!(G_dict, j, v, mkt, sp) < sp.infty &&
        if G_recursor!(G_dict, j, v, mkt, sp) < G_recursor!(G_dict, j - 1, v, mkt, sp)
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

        issorted(mkt.perm) || @warn "t not sorted; table rows will differ from input"
        display(G_table)
    end

    # In the valuation we just use the identity permutation as invp
    # To prevent wasteful permuting and then invpermuting
    return mkt.perm[X], valuation_nopermute(X, mkt)
end
