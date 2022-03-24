using OptimalApplication
using DataFrames
using PrettyTables
using Random
using StatsBase
using Statistics
using Base.Threads
import Printf: @sprintf
using BenchmarkTools
using UnicodePlots

const fullscale = false

# A long benchmark; tweak parameters with caution.
# Set fullscale = false to run a smaller benchmark to check formatting etc.

# Number of times to repeat each computation, where min of these is reported as time
BenchmarkTools.DEFAULT_PARAMETERS.samples = fullscale ? 5 : 2
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1 # ... which is the default

# Number of markets to test at each intersection of the experimental variables
const n_markets = fullscale ? 20 : 2
const bnbcutoff = fullscale ? 33 : 6
const twottbnbcutoff = 2^bnbcutoff

# Sizes of markets to test
const marketsizes_SCM = fullscale ? 4 .^ (2:7) : [5, 10]
const marketsizes_VCM = fullscale ? 2 .^ (3:9) : [5, 10]

function printheader(s)
    printstyled(s * "\n", bold=true, color=222)
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

    sizes = zeros(Int, n_markets, length(marketsizes_SCM))
    times_list = zeros(Float64, n_markets, length(marketsizes_SCM))
    times_heap = zeros(Float64, n_markets, length(marketsizes_SCM))

    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for (j, m) in enumerate(marketsizes_SCM)
            mkt = randSCM(m)
            sizes[i, j] = m
            times_list[i, j] = @belapsed applicationorder_list($mkt)
            times_heap[i, j] = @belapsed applicationorder_heap($mkt)
        end
    end

    df = DataFrame(
        "m" => sizes[:],
        "time_list" => times_list[:],
        "time_heap" => times_heap[:]
    )

    plts = Plot[]
    
    large_idx = df[!, :m] .== marketsizes_SCM[end]
    for c in names(df)[2:end]
        push!(
            plts,
            histogram(
                df[large_idx, c],
                ylabel=c,
                title="Time when m = $(marketsizes_SCM[end])"
            )
        )
    end

    kv_pairs = Pair[]
    for c in names(df)
        if c != "m"
            push!(kv_pairs, c => mean)
            push!(kv_pairs, c => std)
        end
    end
    return combine(groupby(df, :m), kv_pairs...), plts, df
end


function benchmark2()
    printheader("Benchmark 2: Heterogeneous-cost algorithms")
    epsilons = [0.5, 0.05]

    dtype = Union{Float64,Missing}

    sizes = zeros(Int, n_markets, length(marketsizes_VCM))
    times_bnb = zeros(dtype, n_markets, length(marketsizes_VCM))
    times_dp = zeros(dtype, n_markets, length(marketsizes_VCM))
    times_fptas = [zeros(dtype, n_markets, length(marketsizes_VCM)) for k in epsilons]
    fill!.((times_bnb, times_dp), missing)
    fill!.(times_fptas, missing)


    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for (j, m) in enumerate(marketsizes_VCM)
            mkt = randVCM(m)
            sizes[i, j] = m
            if m ≤ bnbcutoff
                times_bnb[i, j] =
                    @belapsed optimalportfolio_branchbound($mkt; maxit=$twottbnbcutoff)
            end
            times_dp[i, j] = @belapsed optimalportfolio_dynamicprogram($mkt)

            for (k, epsilon) in enumerate(epsilons)
                times_fptas[k][i, j] = @belapsed optimalportfolio_fptas($mkt, $epsilon)
            end
        end
    end

    df = DataFrame(
        "m" => sizes[:],
        "time_bnb" => times_bnb[:],
        "time_dp" => times_dp[:],
        ["time_fptas_$epsilon" => times_fptas[k][:] for (k, epsilon) in enumerate(epsilons)]...
    )

    plts = Plot[]

    large_idx = df[!, :m] .== marketsizes_VCM[end]
    for c in setdiff(names(df), ("m", "time_bnb"))
        push!(
            plts,
            histogram(
                Array{Float64}(df[large_idx, c]),
                ylabel=c,
                title="Time when m = $(marketsizes_SCM[end])"
            )
        )
    end

    kv_pairs = Pair[]
    for c in names(df)
        if c != "m"
            push!(kv_pairs, c => mean)
            push!(kv_pairs, c => std)
        end
    end
    return combine(groupby(df, :m), kv_pairs...), plts, df
end


function fmter(v, i, j)
    j == 1 && return v

    v = ismissing(v) ? "—" : @sprintf "%1.2f" 1000 * v

    return isodd(j) ? "($v)" : v
end


@time begin
    println()
    bm1, plts1, df1 = benchmark1()
    display(pretty_table(bm1, formatters=fmter))
    for pl in plts1 display(pl) end
    println("\n")
    bm2, plts2, df2 = benchmark2()
    display(pretty_table(bm2, formatters=fmter))
    for pl in plts2 display(pl) end
    println("\n")
end
