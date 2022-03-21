using OptimalApplication
using DataFrames
using PrettyTables
using Random
using StatsBase
using Statistics
using Base.Threads
import Printf: @sprintf

fullscale = true

# A long benchmark; tweak parameters with caution.
# Set fullscale = false to run a smaller benchmark to check formatting etc.

# Number of times to repeat each computation, where min of these is reported as time
n_reps = fullscale ? 3 : 2

# Number of markets to test at each intersection of the experimental variables
n_markets = fullscale ? 20 : 2

function printheader(s)
    printstyled(s * "\n", bold = true, color = 222)
end

function make_correlated_market(m)
    t = ceil.(Int, 10 * randexp(m))
    sort!(t)
    f = inv.(t .+ 10 * rand(m))
    g = rand(5:10, m)
    H = sum(g) ÷ 2
    return f, t, g, H
end

randVCM(m) = VariedCostsMarket(make_correlated_market(m)...)
function randSCM(m)
    f, t, _, _ = make_correlated_market(m)
    return SameCostsMarket(f, t, m ÷ 2)
end


function benchmark1()
    printheader("Benchmark 1: Homogeneous-cost algorithms")
    M = fullscale ? 4 .^ (2:6) : [5, 10]

    sizes = zeros(Int64, n_markets, length(M))
    times_list = zeros(Float64, n_markets, length(M))
    times_heap = zeros(Float64, n_markets, length(M))

    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for (j, m) in enumerate(M)
            mkt = randSCM(m)
            sizes[i, j] = m
            times_list[i, j] = minimum(@elapsed applicationorder_list(mkt) for r in 1:n_reps)
            times_heap[i, j] = minimum(@elapsed applicationorder_heap(mkt) for r in 1:n_reps)
        end
    end

    df = DataFrame("m" => sizes[:], "time_list" => times_list[:], "time_heap" => times_heap[:])

    kv_pairs = Pair[]
    for c in names(df)
        if c != "m"
            push!(kv_pairs, c => mean)
            push!(kv_pairs, c => std)
        end
    end
    return combine(groupby(df, :m), kv_pairs...), df
end


function benchmark2()
    printheader("Benchmark 2: Heterogeneous-cost algorithms")
    M = fullscale ? 2 .^ (3:9) : [5, 10]
    bnbcutoff = fullscale ? 33 : 6
    epsilons = [0.5, 0.05]

    dtype = Union{Float64,Missing}

    sizes = zeros(Int64, n_markets, length(M))
    times_bnb = zeros(dtype, n_markets, length(M))
    times_dp = zeros(dtype, n_markets, length(M))
    times_fptas = [zeros(dtype, n_markets, length(M)) for k in epsilons]
    fill!.((times_bnb, times_dp), missing)
    fill!.(times_fptas, missing)


    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for (j, m) in enumerate(M)
            mkt = randVCM(m)
            sizes[i, j] = m
            if m ≤ bnbcutoff
                times_bnb[i, j] = minimum(@elapsed optimalportfolio_branchbound(mkt; maxit = 40000) for r in 1:n_reps)
            end
            times_dp[i, j] = minimum(@elapsed optimalportfolio_dynamicprogram(mkt) for r in 1:n_reps)

            for (k, epsilon) in enumerate(epsilons)
                times_fptas[k][i, j] = minimum(@elapsed optimalportfolio_fptas(mkt, epsilon) for r in 1:n_reps)
            end
        end
    end

    df = DataFrame(
        "m" => sizes[:],
        "time_bnb" => times_bnb[:],
        "time_dp" => times_dp[:],
        ["time_fptas_$epsilon" => times_fptas[k][:] for (k, epsilon) in enumerate(epsilons)]...)

    kv_pairs = Pair[]
    for c in names(df)
        if c != "m"
            push!(kv_pairs, c => mean)
            push!(kv_pairs, c => std)
        end
    end
    return combine(groupby(df, :m), kv_pairs...), df
end


function fmter(v, i, j)
    j == 1 && return v

    v = ismissing(v) ? "—" : @sprintf "%1.2f" 1000*v

    return isodd(j) ? "($v)" : v
end


@time begin
    println()
    bm1 = benchmark1()
    display(pretty_table(bm1[1], formatters = fmter))
    println("\n")
    bm2 = benchmark2()
    display(pretty_table(bm2[1], formatters = fmter))
    println("\n")
end
