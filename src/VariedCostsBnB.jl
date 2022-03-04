mutable struct Node{jType<:Integer}
    I::Set{jType}
    N::Set{jType}
    t̄::Dict{jType,<:AbstractFloat}
    H̄::Real
    v_I::AbstractFloat
    v_LP::AbstractFloat
    children::Set{UInt64}
    isleaf::Bool

    function Node(I::Set, N::Set, t̄::Dict, H̄::Real, v_I::Real, mkt::VariedCostsMarket)
        # For generating a new node with its LP relaxation value and empty child set

        # Leaf node
        if isempty(N)
            return new{eltype(I)}(I, N, t̄, H̄, v_I, v_I, Set{UInt64}(), true)
        end

        j_order = collect(N)
        j_order[:] = j_order[sortperm(collect(mkt.f[j] * t̄[j] / mkt.g[j] for j in j_order), rev = true)]

        v_LP = v_I
        H_left = H̄
        for j in j_order
            if mkt.g[j] ≤ H_left
                H_left -= mkt.g[j]
                v_LP += mkt.f[j] * t̄[j]
            else
                v_LP += mkt.f[j] * t̄[j] * H_left / mkt.g[j]
                @goto done
            end
        end
        # Need to use goto here for the case in which all of N fits within H̄
        @label done

        return new{eltype(I)}(I, N, t̄, H̄, v_I, v_LP, Set{UInt64}(), false)
    end
end

hash(nd::Node) = hash((nd.I, nd.N))
isless(nd1::Node, nd2::Node) = isless(nd1.v_LP, nd2.v_LP)



"""
    generatechildren!(nd, mkt)

With respect to the data in `mkt`, generates the child(ren) of node `nd` and writes 
their hashes to `nd.children`.
"""
function generatechildren!(nd::Node, mkt::VariedCostsMarket)
    newI = copy(nd.I)
    newN = copy(nd.N)

    fltr = filter(j -> mkt.g[j] ≤ nd.H̄, nd.N)

    # No way to add any school to this: just skip to the leaf node
    if isempty(fltr)
        child = Node(newI, Set{eltype(nd.I)}(), nd.t̄, nd.H̄, nd.v_I, mkt)
        push!(nd.children, hash(child))
        return (child,)
    end

    # School we will branch on.
    # You may be tempted to divide by g[j] here, but this is unwise, 
    # because the optimal addition on this branch is always the affordable
    # school with highest ft[j]. Different from knapsack.
    i = argmax(j -> mkt.f[j] * nd.t̄[j], fltr)

    delete!(newN, i)

    # Node without i
    t̄2 = copy(nd.t̄)
    delete!(t̄2, i)

    child2 = Node(newI, newN, t̄2, nd.H̄, nd.v_I, mkt)

    if mkt.g[i] < nd.H̄
        # Node with i
        t̄1 = copy(nd.t̄)
        for j in keys(t̄1)
            if t̄1[j] ≤ t̄1[i]
                t̄1[j] *= mkt.omf[i]
            else
                t̄1[j] -= mkt.f[i] * nd.t̄[i]
            end
        end
        delete!(t̄1, i)

        child1 = Node(union(newI, i), newN, t̄1, nd.H̄ - mkt.g[i],
            nd.v_I + mkt.f[i] * nd.t̄[i],
            mkt)

        push!(nd.children, hash(child1), hash(child2))
        return child1, child2
    else
        # Then mkt.g[i] == nd.H̄
        # Leaf node with i
        t̄1 = copy(nd.t̄)
        for j in keys(t̄1)
            if t̄1[j] ≤ t̄1[i]
                t̄1[j] *= mkt.omf[i]
            else
                t̄1[j] -= mkt.f[i] * nd.t̄[i]
            end
        end
        delete!(t̄1, i)

        child1 = Node(union(newI, i), Set{eltype(nd.I)}(), t̄1, 0, nd.v_I + mkt.f[i] * nd.t̄[i], mkt)

        push!(nd.children, hash(child1), hash(child2))
        return child1, child2
    end
end


function fathom!(tree, fathomhash, LB_hash)
    for ndhash_ in tree[fathomhash].children
        # Need to check if it has a key since we might have fathomed it somewhere else
        haskey(tree, ndhash_) && fathom!(tree, ndhash_, LB_hash) # && !tree[ndhash_].isleaf
    end

    delete!(tree, fathomhash)
end


function optimalportfolio_branchbound(mkt; maxit = 10000::Integer, verbose = false::Bool)
    mkt.m > 32 && @warn "Branch and bound is slow for large markets"

    tree = Dict{UInt64,Node}()

    C = Set(1:length(mkt.t))

    rootnode = Node(Set{eltype(C)}(), C, Dict(zip(1:mkt.m, Float64.(mkt.t))), mkt.H, 0, mkt)
    LB = 0
    LB_hash = hash(rootnode)
    push!(tree, LB_hash => rootnode)

    for k in 1:maxit
        verbose && @show k, LB, length(tree)

        fltr = filter(hash_nd -> !hash_nd[2].isleaf && isempty(hash_nd[2].children), tree)

        if isempty(fltr)
            # All branches have either reached leaves or fathomed: done
            return collect(tree[LB_hash].I), tree[LB_hash].v_I
        else
            # Select the node with the highest UB
            thisnodehash = argmax(fltr)
        end

        children = generatechildren!(tree[thisnodehash], mkt)

        for child in children
            childhash = hash(child)

            push!(tree, childhash => child)
            if child.v_I > LB
                LB = child.v_I
                LB_hash = childhash
            end

            fathomhash = findfirst(nd -> nd.v_LP < LB, tree)
            while !isnothing(fathomhash)
                verbose && println("Fathoming node $fathomhash")
                fathom!(tree, fathomhash, LB_hash)

                fathomhash = findfirst(nd -> nd.v_LP < LB, tree)
            end
        end
    end

    @warn "Unable to find optimum in $maxit iterations; returning best so far.\n" *
        "         Worst-case optimality ratio: $(tree[LB_hash].v_I/rootnode.v_LP)"
    return collect(tree[LB_hash].I), tree[LB_hash].v_I
end