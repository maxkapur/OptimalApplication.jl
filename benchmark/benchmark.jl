using OptimalApplication
using DataFrames
using PrettyTables
using Random
using StatsBase
using Statistics
using Base.Threads
import Printf: @sprintf
# using BenchmarkTools
using UnicodePlots

const fullscale = true

# A long benchmark; tweak parameters with caution.
# Set fullscale = false to run a smaller benchmark to check formatting etc.

# Number of times to repeat each computation, where min of these is reported as time
const n_reps = fullscale ? 5 : 2
# BenchmarkTools.DEFAULT_PARAMETERS.samples = n_reps
# BenchmarkTools.DEFAULT_PARAMETERS.evals = 1 # ... which is the default

# Number of markets to test at each intersection of the experimental variables
const n_markets = fullscale ? 40 : 2
const bnbcutoff = fullscale ? 33 : 6
const twottbnbcutoff = 2^bnbcutoff

# Sizes of markets to test
const marketsizes_SCM = fullscale ? 4 .^ (2:7) : [5, 10]
const marketsizes_VCM = fullscale ? 2 .^ (3:9) : [5, 10]

function makecorrelatedmarketdata(m)
    t = ceil.(Int, 10 * randexp(m))
    sort!(t)
    f = inv.(t .+ 10 * rand(m))
    g = rand(5:10, m)
    H = sum(g) ÷ 2
    return f, t, g, H
end

randVCM(m) = VariedCostsMarket(makecorrelatedmarketdata(m)...)

function randSCM(m)
    f, t, _, _ = makecorrelatedmarketdata(m)
    return SameCostsMarket(f, t, m ÷ 2)
end

const mkts_SCM = SameCostsMarket[randSCM(m) for i in 1:n_markets, m in marketsizes_SCM]
const mkts_VCM = VariedCostsMarket[randVCM(m) for i in 1:n_markets, m in marketsizes_VCM]


function makeunicodeplots(df::DataFrame)
    large_idx = df[!, :m] .== maximum(df[!, :m])
    return [
        histogram(
            Vector{Float64}(df[large_idx, c]),
            ylabel=c,
            title="Time when m = $(marketsizes_SCM[end])"
        )
        for c in setdiff(names(df), ("m", "time_bnb"))
    ]
end


function collectmeanstd(df::DataFrame)
    kv_pairs = Pair[]
    for c in names(df)
        if c != "m"
            push!(kv_pairs, c => mean)
            push!(kv_pairs, c => std)
        end
    end
    return combine(groupby(df, :m), kv_pairs...)
end


function printheader(s)
    printstyled(s * "\n", bold=true, color=222)
end


function benchmark1()
    printheader("Benchmark 1: Homogeneous-cost algorithms")

    sizes = Int[m for _ in 1:n_markets, m in marketsizes_SCM]
    times_list = fill(Inf, n_markets, length(marketsizes_SCM))
    times_heap = fill(Inf, n_markets, length(marketsizes_SCM))

    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for j in 1:length(marketsizes_SCM), _ in 1:n_reps
            times_list[i, j] = min(times_list[i, j], @elapsed applicationorder_list(mkts_SCM[i, j]))
            times_heap[i, j] = min(times_heap[i, j], @elapsed applicationorder_heap(mkts_SCM[i, j]))
        
            # times_list[i, j] = @belapsed applicationorder_list($(mkts_SCM[i, j]))
            # times_heap[i, j] = @belapsed applicationorder_heap($(mkts_SCM[i, j]))
        end
    end

    df = DataFrame(
        "m" => sizes[:],
        "time_list" => times_list[:],
        "time_heap" => times_heap[:]
    )

    plts = makeunicodeplots(df)
    meanstd = collectmeanstd(df)

    return meanstd, plts, df
end


function benchmark2()
    printheader("Benchmark 2: Heterogeneous-cost algorithms")
    epsilons = [0.5, 0.05]

    sizes = Int[m for _ in 1:n_markets, m in marketsizes_VCM]
    times_bnb = fill(Inf, n_markets, length(marketsizes_VCM))
    times_dp = fill(Inf, n_markets, length(marketsizes_VCM))
    times_fptas = fill(Inf, n_markets, length(marketsizes_VCM), length(epsilons))

    for i in 1:n_markets
        println("  i = $i of $n_markets")
        for (j, m) in enumerate(marketsizes_VCM), _ in 1:n_reps
            if m ≤ bnbcutoff
                times_bnb[i, j] = min(times_bnb[i, j], @elapsed optimalportfolio_branchbound(mkts_VCM[i, j]; maxit=twottbnbcutoff))
                # times_bnb[i, j] =
                #     @belapsed optimalportfolio_branchbound($(mkts_VCM[i, j]); maxit=$twottbnbcutoff)
            end

            times_dp[i, j] = min(times_dp[i, j], @elapsed optimalportfolio_dynamicprogram(mkts_VCM[i, j]))
            # times_dp[i, j] = @belapsed optimalportfolio_dynamicprogram($(mkts_VCM[i, j]))

            for (k, epsilon) in enumerate(epsilons)
                times_fptas[i, j, k] = min(times_fptas[i, j, k], @elapsed optimalportfolio_fptas(mkts_VCM[i, j], epsilon))
                # times_fptas[i, j, k] = @belapsed optimalportfolio_fptas($(mkts_VCM[i, j]), $epsilon)
            end
        end
    end

    df = DataFrame(
        "m" => sizes[:],
        "time_bnb" => times_bnb[:],
        "time_dp" => times_dp[:],
        ["time_fptas_$epsilon" => times_fptas[:, :, k][:] for (k, epsilon) in enumerate(epsilons)]...
    )

    plts = makeunicodeplots(df)
    meanstd = collectmeanstd(df)

    return meanstd, plts, df
end


function fmter(v, i, j)
    j == 1 && return v

    v = ismissing(v) || isinf(v) || isnan(v) ? "—" : @sprintf "%1.2f" 1000 * v

    return isodd(j) ? "($v)" : v
end


function doeverything()
    println()
    bm1, plts1, = benchmark1()
    display(pretty_table(bm1, formatters=fmter))
    for pl in plts1
        display(pl)
        println()
    end
    bm2, plts2, = benchmark2()
    display(pretty_table(bm2, formatters=fmter))
    for pl in plts2
        display(pl)
        println()
    end
end


@time doeverything()
