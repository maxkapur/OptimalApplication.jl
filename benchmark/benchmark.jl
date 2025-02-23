using OptimalApplication
using DataFrames
using PrettyTables
using Random
using StatsBase
using Statistics
using Base.Threads
import Printf: @sprintf
using UnicodePlots

const fullscale = false

# A long benchmark; tweak parameters with caution.
# Set fullscale = false to run a smaller benchmark to check formatting etc.

# Number of times to repeat each computation, where min of these is reported as time
const n_reps = fullscale ? 3 : 2

# Number of markets to test at each intersection of the experimental variables
const n_markets = fullscale ? 50 : 5

# Cutoffs to exclude certain large computations
const bnbcutoff = fullscale ? 33 : 10
const slowdpcutoff = fullscale ? 33 : 10
const twottbnbcutoff = 2^bnbcutoff
const fptascutoff_m = fullscale ? 300 : 200
const fptascutoff_eps = fullscale ? 0.0 : 0.1

# Market sizes
marketsizes_SCM = fullscale ? 4 .^ (2:7) : 2 .^ (3:6)
marketsizes_VCM = fullscale ? 2 .^ (3:8) : 2 .^ (3:6)

function makeunicodeplots(df::DataFrame)
    large_idx = df[!, :m] .== maximum(df[!, :m])
    return [
        histogram(
            Vector{Float64}(df[large_idx, c]),
            ylabel = c,
            title = "Time when m = $(marketsizes_SCM[end])",
            canvas = BlockCanvas,
        ) for c in setdiff(names(df), ("m", "time_bnb")) if !(Inf in df[!, c])
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


function benchmark1(marketsizes_SCM = marketsizes_SCM)
    @info "Benchmark 1: Homogeneous-cost algorithms"
    mkts_SCM = SameCostsMarket[SameCostsMarket(m) for m in marketsizes_SCM, i = 1:n_markets]

    sizes = Int[m for m in marketsizes_SCM, _ = 1:n_markets]
    times_list = fill(Inf, length(marketsizes_SCM), n_markets)
    times_heap = fill(Inf, length(marketsizes_SCM), n_markets)

    @threads for j = 1:n_markets
        println("  j = ", @sprintf("%2d", j), " of $n_markets")
        for (i, _) in enumerate(marketsizes_SCM), _ = 1:n_reps
            times_list[i, j] =
                min(times_list[i, j], @elapsed applicationorder_list(mkts_SCM[i, j]))
            times_heap[i, j] =
                min(times_heap[i, j], @elapsed applicationorder_heap(mkts_SCM[i, j]))
        end
    end

    df = DataFrame(
        "m" => sizes[:],
        "time_list" => times_list[:],
        "time_heap" => times_heap[:],
    )

    plts = makeunicodeplots(df)
    meanstd = collectmeanstd(df)

    return meanstd, plts, df
end


function benchmark2(marketsizes_VCM = marketsizes_VCM)
    @info "Benchmark 2: Heterogeneous-cost algorithms"
    delta = 0.1             # for slow DP
    epsilons = [0.5, 0.05]  # for FPTAS
    mkts_VCM =
        VariedCostsMarket[VariedCostsMarket(m) for m in marketsizes_VCM, i = 1:n_markets]

    sizes = Int[m for m in marketsizes_VCM, _ = 1:n_markets]
    times_bnb = fill(Inf, length(marketsizes_VCM), n_markets)
    times_dp = fill(Inf, length(marketsizes_VCM), n_markets)
    times_dp_slow = fill(Inf, length(marketsizes_VCM), n_markets)
    times_fptas = fill(Inf, length(epsilons), length(marketsizes_VCM), n_markets)

    @threads for j = 1:n_markets
        println("  j = ", @sprintf("%2d", j), " of $n_markets")
        for (i, m) in enumerate(marketsizes_VCM), _ = 1:n_reps
            if m ≤ bnbcutoff
                times_bnb[i, j] = min(
                    times_bnb[i, j],
                    @elapsed optimalportfolio_branchbound(
                        mkts_VCM[i, j];
                        maxit = twottbnbcutoff,
                    )
                )
                # times_bnb[i, j] =
                #     @belapsed optimalportfolio_branchbound($(mkts_VCM[i, j]); maxit=$twottbnbcutoff)
            end

            times_dp[i, j] = min(
                times_dp[i, j],
                @elapsed optimalportfolio_dynamicprogram(mkts_VCM[i, j])
            )
            # times_dp[i, j] = @belapsed optimalportfolio_dynamicprogram($(mkts_VCM[i, j]))

            if m ≤ slowdpcutoff
                times_dp_slow[i, j] = min(
                    times_dp_slow[i, j],
                    @elapsed optimalportfolio_dynamicprogram_slow(mkts_VCM[i, j], delta)
                )
            end

            for (k, epsilon) in enumerate(epsilons)
                if m ≤ fptascutoff_m || epsilon > fptascutoff_eps
                    times_fptas[k, i, j] = min(
                        times_fptas[k, i, j],
                        @elapsed optimalportfolio_fptas(mkts_VCM[i, j], epsilon)
                    )
                    # times_fptas[i, j, k] = @belapsed optimalportfolio_fptas($(mkts_VCM[i, j]), $epsilon)
                end
            end
        end
    end

    df = DataFrame(
        "m" => sizes[:],
        "time_bnb" => times_bnb[:],
        "time_dp" => times_dp[:],
        "time_dp_slow_$delta" => times_dp_slow[:],
        [
            "time_fptas_$epsilon" => times_fptas[k, :, :][:] for
            (k, epsilon) in enumerate(epsilons)
        ]...,
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
    @info "Benchmark parameters" fullscale n_markets n_reps nthreads()
    println()
    bm1, plts1, = benchmark1()
    display(pretty_table(bm1, formatters = fmter))
    for pl in plts1
        display(pl)
        println()
    end
    bm2, plts2, = benchmark2()
    display(pretty_table(bm2, formatters = fmter))
    for pl in plts2
        display(pl)
        println()
    end
end


@time doeverything()
