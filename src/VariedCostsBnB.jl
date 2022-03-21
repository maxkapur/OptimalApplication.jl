"""
Contains information about a subproblem in the branch and bound scheme.
`I` is a set of schools that are "in" the portfolio and `N` is the set that 
are "negotiable" and used to generate an LP upper bound and child nodes.
"""
mutable struct Node
    I::Set{Int}
    N::Set{Int}
    t̄::Dict{Int,Float64}
    H̄::Int
    v_I::Float64
    v_LP::Float64

    function Node(I::Set, N::Set, t̄::Dict, H̄::Int, v_I::Float64, mkt::VariedCostsMarket)
        # For generating a new node with its LP relaxation value and empty child set

        # Leaf node
        if isempty(N)
            return new(I, N, t̄, H̄, v_I, v_I)
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

        return new(I, N, t̄, H̄, v_I, v_LP)
    end
end

hash(nd::Node) = hash((nd.I, nd.N))
# Depth-first strategy:
# isless(nd1::Node, nd2::Node) = isless(nd1.v_I, nd2.v_I)
# Breadth-first strategy:
# isless(nd1::Node, nd2::Node) = isless(nd1.v_LP, nd2.v_LP)



"""
    generatechildren(nd, mkt)

With respect to the data in `mkt`, generates the child(ren) of node `nd` and writes 
their hashes to `nd.children`.
"""
function generatechildren(nd::Node, mkt::VariedCostsMarket)
    fltr = filter(j -> mkt.g[j] ≤ nd.H̄, nd.N)

    # No way to add any school to this: just skip to the leaf node
    if isempty(fltr)
        child = Node(nd.I, Set{Int}(), nd.t̄, nd.H̄, nd.v_I, mkt)
        return (child,)
    end

    # School we will branch on. In principle it can by any school in fltr.
    # i = argmax(j -> mkt.f[j] * nd.t̄[j] / mkt.g[j], fltr)
    i = argmax(j -> mkt.f[j] * nd.t̄[j], fltr)

    newN = setdiff(nd.N, i)

    # Node without i
    t̄2 = copy(nd.t̄)
    delete!(t̄2, i)

    child2 = Node(nd.I, newN, t̄2, nd.H̄, nd.v_I, mkt)

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

        child1 = Node(union(nd.I, i), newN, t̄1, nd.H̄ - mkt.g[i],
            nd.v_I + mkt.f[i] * nd.t̄[i],
            mkt)

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

        child1 = Node(union(nd.I, i), Set{Int}(), t̄1, 0, nd.v_I + mkt.f[i] * nd.t̄[i], mkt)

        return child1, child2
    end
end


"""
    optimalportfolio_branchbound(mkt::VariedCostsMarket; maxit=10000, verbose=false)

Use the branch-and-bound algorithm to produce the optimal portfolio for the
market `mkt` with varying application costs. Intractable for large markets. 
"""
function optimalportfolio_branchbound(mkt; maxit = 100000::Int, verbose = false::Bool)
    mkt.m ≥ 33 && @warn "Branch and bound is slow for large markets"

    C = Set(1:length(mkt.t))

    rootnode = Node(Set{Int}(), C, Dict(zip(1:mkt.m, Float64.(mkt.t))), mkt.H, 0.0, mkt)
    LB = 0
    LB_node = rootnode

    # Commented out is an implementation that stores tree as a heap. The dict is faster because
    # fathoming nodes is slow on the heap. 

    # tree = MutableBinaryMaxHeap{Node{eltype(C)}}()        
    # treekeys = Set{Int64}()

    # push!(treekeys, push!(tree, rootnode))

    tree = Dict{UInt,Node}()
    push!(tree, hash(rootnode) => rootnode)

    for k in 1:maxit
        verbose && @show k, LB, length(tree)

        if isempty(tree)
            # All branches have either reached leaves or fathomed: done
            return collect(LB_node.I), LB
        else
            # Select the node with the highest UB
            # thisnode, thisnodehandle = top_with_handle(tree)
            # pop!(tree)
            # delete!(treekeys, thisnodehandle)
        
            thisnodehandle, thisnode = argmax(hash_nd -> hash_nd[2].v_LP, tree)
            delete!(tree, thisnodehandle)
            
            # Another option: Select the node with best obj value. Works pretty bad. 
            # thisnodehandle, thisnode = argmax(hash_nd -> hash_nd[2].v_I, tree)
            # delete!(tree, thisnodehandle)
        end

        children = generatechildren(thisnode, mkt)

        newLB = false
        for child in children
            if child.v_I > LB
                newLB = true
                LB_node = child
                LB = child.v_I
            end

            # isempty(child.N) || push!(treekeys, push!(tree, child))
            if child.v_LP > LB && !isempty(child.N)
                push!(tree, hash(child) => child)
            end
        end

        if newLB
            # for k in treekeys
            for k in keys(tree)
                if tree[k].v_LP < LB
                    verbose && println("Fathoming node $k")
                    delete!(tree, k)
                    # delete!(treekeys, k)
                end
            end
        end
    end

    @warn "Unable to find optimum in $maxit iterations; returning best so far.\n" *
          "         Worst-case optimality ratio: $(LB/rootnode.v_LP)"
    return collect(LB_node.I), LB
end
