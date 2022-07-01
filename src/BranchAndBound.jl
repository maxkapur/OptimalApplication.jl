"""
Contains information about a subproblem in the branch and bound scheme.
`I` is a set of schools that are "in" the portfolio and `N` is the set that 
are "negotiable" and used to generate an LP upper bound and child nodes.
"""
mutable struct Node
    I::Vector{Int}
    N::Vector{Int}
    t_bar::Vector{Float64}
    H_bar::Int
    v_I::Float64              # Valuation of the schools in I
    v_LP::Float64             # LP upper bound

    "Generate a new node with its LP relaxation value and empty child set."
    function Node(
        I::Vector{Int},
        N::Vector{Int},
        t_bar::Vector{Float64},
        H_bar::Int,
        v_I::Float64,
        mkt::VariedCostsMarket
    )
        # Leaf node
        isempty(N) && new(I, N, t_bar, H_bar, v_I, v_I)

        sort!(
            N,
            by=function (j)
                mkt.f[j] * t_bar[j] / mkt.g[j]
            end,
            rev=true,
            alg=InsertionSort # Typically just one item is out of place
        )

        v_LP = v_I
        H_left = H_bar
        @inbounds for j in N
            if mkt.g[j] ≤ H_left
                H_left -= mkt.g[j]
                v_LP += mkt.f[j] * t_bar[j]
            else
                v_LP += mkt.f[j] * t_bar[j] * H_left / mkt.g[j]
                break
            end
        end

        return new(I, N, t_bar, H_bar, v_I, v_LP)
    end
end


isless(n::Node, m::Node) = isless(n.v_LP, m.v_LP)


"""
    generatechildren(nd, mkt)

With respect to the data in `mkt`, generates the child(ren) of node `nd`.
"""
function generatechildren(nd::Node, mkt::VariedCostsMarket)::Vector{Node}
    newN = filter(j -> mkt.g[j] ≤ nd.H_bar, nd.N)

    # No way to add any school to this: just skip to the leaf node
    if isempty(newN)
        child = Node(nd.I, Int[], nd.t_bar, nd.H_bar, nd.v_I, mkt)
        return [child]
    end

    # School we will branch on. In principle it can by any school in newN.
    i = popfirst!(newN) # = argmax(j -> mkt.f[j] * nd.t_bar[j], newN) by sortation policy

    # Node without i
    child2 = Node(nd.I, newN, nd.t_bar, nd.H_bar, nd.v_I, mkt)

    if mkt.g[i] < nd.H_bar
        # Node with i
        t_bar1 = copy(nd.t_bar)
        @inbounds for j in newN
            if t_bar1[j] ≤ nd.t_bar[i]
                t_bar1[j] *= (1 - mkt.f[i])
            else
                t_bar1[j] -= mkt.f[i] * nd.t_bar[i]
            end
        end

        child1 = Node(
            vcat(nd.I, i),
            newN,
            t_bar1,
            nd.H_bar - mkt.g[i],
            nd.v_I + mkt.f[i] * nd.t_bar[i],
            mkt
        )

        return [child1, child2]
    else
        child1 = Node(
            vcat(nd.I, i),
            Int[],
            Float64[], # nd.t_bar
            0,
            nd.v_I + mkt.f[i] * nd.t_bar[i],
            mkt
        )

        return [child1, child2]
    end
end


"""
    optimalportfolio_branchbound(mkt::VariedCostsMarket; maxit=10000, verbose=false) -> X, v

Use the branch-and-bound algorithm to produce the optimal portfolio for the `VariedCostsMarket`
defined by `mkt`. Intractable for large markets. 

```julia-repl
julia> mkt = VariedCostsMarket([0.2, 0.5, 0.1, 0.6, 0.1], [1, 4, 9, 1, 8], [2, 4, 2, 5, 1], 8);

julia> optimalportfolio_branchbound(mkt)
([2, 3, 5], 3.24)
```
"""
function optimalportfolio_branchbound(mkt::VariedCostsMarket; maxit::Int=100000, verbose::Bool=false)::Tuple{Vector{Int},Float64}
    mkt.m ≥ 33 && @warn "Branch and bound is slow for large markets"

    C = collect(1:mkt.m)

    rootnode::Node = Node(Int[], C, Float64.(mkt.t), mkt.H, 0.0, mkt)
    LB::Float64 = 0.0
    LB_node::Node = rootnode
    tree = BinaryMaxHeap{Node}([rootnode])

    for k in 1:maxit
        verbose && @show k, LB, length(tree)

        # All branches have either reached leaves or fathomed: done
        isempty(tree) && return mkt.perm[collect(LB_node.I)], LB

        # Inspect the node with the highest UB
        children::Vector{Node} = generatechildren(pop!(tree), mkt)

        for child in children
            if child.v_LP > LB
                if child.v_I > LB
                    LB_node = child
                    LB = child.v_I
                end
                isempty(child.N) || push!(tree, child)
            end
        end
    end

    @warn "Unable to find optimum in $maxit iterations; returning best so far.\n" *
          "         Worst-case optimality ratio: $(LB/rootnode.v_LP)"
    return mkt.perm[collect(LB_node.I)], LB
end
