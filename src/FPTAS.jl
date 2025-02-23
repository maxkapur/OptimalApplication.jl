# Used by fptas below
# TODO: Consider copying mkt.f here for memory localization
struct ScaleParams
    t::Vector{Int}
    # ft::Vector{Float64} = mkt.f .* t
    v_UB::Int
    infty::Int

    function ScaleParams(mkt::VariedCostsMarket, ε::Float64)
        @assert 0 < ε < 1
        P = max(0, ceil(Int, log2(Int(mkt.m)^2 / (ε * sum(mkt.f .* mkt.t)))))
        t = mkt.t .* 2^P
        v_UB = ceil(Int, sum(mkt.f .* t))
        infty = mkt.H + 1

        return new(t, v_UB, infty)
    end
end

# Used by fptas below
@inbounds function G_recursor!(
    G_dict::Dict{Tuple{Int,Int},Int},
    j::Int,
    v::Int,
    mkt::VariedCostsMarket,
    sp::ScaleParams,
)::Int
    return get(G_dict, (j, v)) do
        if v ≤ 0
            return 0
        elseif j == 0 || sp.t[j] < v || v ≥ sp.v_UB
            return sp.infty
        else
            if mkt.f[j] < 1
                # Clamping prevents over/underflow: for any v<0 or v≥v_UB the function
                # is trivially defined, so recording any more extreme number is meaningless
                v_minus_Δ =
                    clamp(ceil(Int, (v - mkt.f[j] * sp.t[j]) / (1 - mkt.f[j])), -1, sp.v_UB)

                push!(
                    G_dict,
                    (j, v) => min(
                        G_recursor!(G_dict, j - 1, v, mkt, sp),
                        mkt.g[j] + G_recursor!(G_dict, j - 1, v_minus_Δ, mkt, sp),
                    ),
                )
                return G_dict[(j, v)]
            else
                push!(
                    G_dict,
                    (j, v) => min(G_recursor!(G_dict, j - 1, v, mkt, sp), mkt.g[j]),
                )
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

```julia-repl
julia> mkt = VariedCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], [2, 4, 2, 5, 1], 8);

julia> optimalportfolio_fptas(mkt, 0.2)
([3, 5, 2], 3.24)
```
"""
function optimalportfolio_fptas(
    mkt::VariedCostsMarket,
    ε::Float64;
    verbose::Bool = false,
)::Tuple{Vector{Int},Float64}
    sp = ScaleParams(mkt, ε)

    G_dict = Dict{Tuple{Int,Int},Int}()
    sizehint!(G_dict, mkt.m * mkt.m ÷ 2)

    # Binary search for max{ w :  G_recursor!(G_dict, mkt.m, w, mkt, sp) ≤ H }
    # In future Julia, may be able to do something like this:

    # vs = 0:sp.v_UB
    # searchsortedlast(vs, by = ...)

    v = 0
    v_current_UB = sp.v_UB
    @inbounds while v + 1 < v_current_UB
        # mid = (v + v_UB) ÷ 2
        mid = midpoint(v, v_current_UB)

        if G_recursor!(G_dict, mkt.m, mid, mkt, sp) > mkt.H
            v_current_UB = mid
        else
            v = mid
        end
    end

    X = Int[]
    @inbounds for j in reverse(1:mkt.m)
        # G_recursor!(G_dict, j, v, mkt, sp) < sp.infty &&
        if G_recursor!(G_dict, j, v, mkt, sp) < G_recursor!(G_dict, j - 1, v, mkt, sp)
            push!(X, j)
            v = clamp(ceil(Int, (v - mkt.f[j] * sp.t[j]) / (1 - mkt.f[j])), -1, sp.v_UB)
        end
    end

    if verbose
        G_table::Array{Int} = fill(sp.infty, (mkt.m, sp.v_UB))
        for (j, v) in keys(G_dict)
            if 0 < j ≤ mkt.m && 0 < v ≤ sp.v_UB
                G_table[j, v] = G_dict[(j, v)]
            end
        end

        issorted(mkt.perm) || @warn "t not sorted; table rows will differ from input"
        display(G_table)
    end

    return mkt.perm[X], valuation_nopermute_sorted(reverse(X), mkt)
end
