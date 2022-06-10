"""
    College(j, f, t, ft, omf)

Contains a college's index `j`, admissions probability `f`, utility value `t`,
their product `ft`, and `1 - f = omf`. Used only by `applicationorder()`.
"""
struct College
    j::Int
    f::Float64
    t::Float64
end

# Overload isless() so that the heap is ordered by expected utility.
isless(c1::College, c2::College) = isless(c1.f * c1.t, c2.f * c2.t)


"""
    applicationorder_list(mkt::SameCostsMarket) -> X, V

Produce the optimal application order `X` and associated valuations `V` 
for the [`SameCostsMarket`](@ref) defined by `mkt`. Uses a list data structure;
typically faster than the equivalent [`applicationorder_heap`](@ref).

All `SameCostsMarket`s satisfy a nestedness property, meaning that the optimal
portfolio when `mkt.h = h` is given by the first `h` entries of the optimal portfolio
when `mkt.h = mkt.m`. For example:

```julia-repl
julia> mkt = SameCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], 5);

julia> x, v = applicationorder_list(mkt)
([2, 3, 5, 4, 1], [2.0, 2.7, 3.24, 3.483, 3.5154])

julia> x[1:4], v[4] # optimal portfolio and valuation for h = 4 
([2, 3, 5, 4], 3.483)
```
"""
function applicationorder_list(
    mkt::SameCostsMarket;
    verbose::Bool=false
)::Tuple{Vector{Int},Vector{Float64}}
    apporder = zeros(Int, mkt.h)
    v = zeros(mkt.h)
    mkt_list = College[College(j, mkt.f[j], mkt.t[j]) for j in 1:mkt.m]

    c_best::College, idx_best::Int = findmax(mkt_list)
    @inbounds for j in 1:mkt.h
        if verbose
            println("Iteration $j")
            println("  ft: ", [mkt_list[i].f * mkt_list[i].t for i in 1:length(mkt_list)])
            println("  Add school $(c_best.j)")
        end

        v[j] = get(v, j - 1, 0) + c_best.f * c_best.t
        apporder[j] = c_best.j

        next_c_best::College = College(0, 1.0, -1.0)
        next_idx_best::Int = 0

        for i in 1:idx_best-1
            mkt_list[i] =
                College(
                    mkt_list[i].j,
                    mkt_list[i].f,
                    mkt_list[i].t * (1 - c_best.f),
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
                    mkt_list[i].t - c_best.f * c_best.t,
                )

            if isless(next_c_best, mkt_list[i-1])
                next_c_best = mkt_list[i-1]
                next_idx_best = i - 1
            end
        end

        c_best, idx_best = next_c_best, next_idx_best

        pop!(mkt_list)
    end

    return mkt.perm[apporder], v
end


"""
    applicationorder_list(mkt::SameCostsMarket) -> X, V

Produce the optimal application order `X` and associated valuations `V` 
for the [`SameCostsMarket`](@ref) defined by `mkt`. Uses a heap data structure;
typically the equivalent [`applicationorder_list`](@ref) is faster.

All `SameCostsMarket`s satisfy a nestedness property, meaning that the optimal
portfolio when `mkt.h = h` is given by the first `h` entries of the optimal portfolio
when `mkt.h = mkt.m`. For example:

```julia-repl
julia> mkt = SameCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], 5);

julia> x, v = applicationorder_heap(mkt)
([2, 3, 5, 4, 1], [2.0, 2.7, 3.24, 3.483, 3.5154])

julia> x[1:4], v[4] # optimal portfolio and valuation for h = 4 
([2, 3, 5, 4], 3.483)
```
"""
function applicationorder_heap(mkt::SameCostsMarket)::Tuple{Vector{Int},Vector{Float64}}
    apporder = zeros(Int, mkt.h)
    v = zeros(mkt.h)

    mkt_heap = BinaryMaxHeap{College}(collect(College(j, mkt.f[j], mkt.t[j]) for j in 1:mkt.m))

    @inbounds for j in 1:mkt.h
        c_k = first(mkt_heap)
        v[j] = get(v, j - 1, 0) + c_k.f * c_k.t
        apporder[j] = c_k.j

        # These two implementations perform almost identically
        mkt_heap = BinaryMaxHeap{College}(
            # College{T}[
            #     College(
            #         c.j,
            #         c.f,
            #         c.j < c_k.j ? c.t * (1 - c_k.f) : c.t - c_k.f * c_k.t,
            #     )
            #     for c in mkt_heap.valtree if c.j != c_k.j
            # ]
            map(filter(c -> c.j != c_k.j, mkt_heap.valtree)) do c
                College(
                    c.j,
                    c.f,
                    c.j < c_k.j ? c.t * (1 - c_k.f) : c.t - c_k.f * c_k.t,
                )
            end
        )
    end

    return mkt.perm[apporder], v
end


# Used by dynamic program below
@inbounds function V_recursor!(
    V_dict::Dict{Tuple{Int,Int},Float64},
    j::Int,
    h::Int,
    mkt::VariedCostsMarket
)::Float64
    return get(V_dict, (j, h)) do
        if j == 0 || h == 0
            return 0.0
        elseif h < mkt.g[j]
            push!(V_dict, (j, h) => V_recursor!(V_dict, j - 1, h, mkt))
            return V_dict[(j, h)]
        else
            push!(V_dict, (j, h) => max(
                V_recursor!(V_dict, j - 1, h, mkt),
                (1 - mkt.f[j]) * V_recursor!(V_dict, j - 1, h - mkt.g[j], mkt) + mkt.f[j] * mkt.t[j]
            ))
            return V_dict[(j, h)]
        end
    end
end


"""
    optimalportfolio_dynamicprogram(mkt::VariedCostsMarket) -> X, v

Use the dynamic program on application costs to produce the optimal portfolio `X` and associated
value `v` for the [`VariedCostsMarket`](@ref) defined by `mkt`. 

```julia-repl
julia> mkt = VariedCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], [2, 4, 2, 5, 1], 8);

julia> optimalportfolio_dynamicprogram(mkt)
([3, 5, 2], 3.24)
```
"""
function optimalportfolio_dynamicprogram(
    mkt::VariedCostsMarket;
    verbose::Bool=false
)::Tuple{Vector{Int},Float64}
    V_dict = Dict{Tuple{Int,Int},Float64}()
    sizehint!(V_dict, mkt.m * mkt.m ÷ 2)

    v = V_recursor!(V_dict, mkt.m, mkt.H, mkt)

    h = mkt.H
    X = Int[]

    for j in reverse(1:mkt.m)
        if V_recursor!(V_dict, j - 1, h, mkt) < V_recursor!(V_dict, j, h, mkt)
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

    return mkt.perm[X], v
end
