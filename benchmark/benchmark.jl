using OptimalApplication
using DataFrames
using Random
using Statistics
using Base.Threads

fullscale = true

# An all-day benchmark; tweak parameters with caution.
# Set fullscale = false to run a smaller benchmark to check formatting etc.

# Number of times to repeat each computation, where min of these is reported as time
n_reps = fullscale ? 3 : 1

function printheader(s)
    printstyled(s*"\n", bold=true, color=222)
end

function make_correlated_market(m)
    A = 10
    t = ceil.(Int, 10 * randexp(m))
    sort!(t)
    f = inv.(t .+ 10 * rand(m))
    g = rand(5:10, m)
    H = sum(g) ÷ 2
    return f, t, g, H
end

function samecosts(ds, M, n_markets)
    sizes = zeros(Int64, n_markets, length(M))
    times = zeros(Float64, n_markets, length(M))

    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for (j, m) in enumerate(M)
            f, t, = make_correlated_market(m)
            sizes[i, j] = m
            times[i, j] = 1000 * minimum(@elapsed applicationorder(f, t; datastructure = ds) for r in 1:n_reps)
        end
    end

    return DataFrame("m" => sizes[:], "time" => times[:])
end


function dp_H(M, n_markets)
    sizes = zeros(Int64, n_markets, length(M))
    times = zeros(Float64, n_markets, length(M))

    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for (j, m) in enumerate(M)
            f, t, g, H = make_correlated_market(m)
            sizes[i, j] = m
            times[i, j] = 1000 * minimum(@elapsed optimalportfolio_dynamicprogram(f, t, g, H) for r in 1:n_reps)
        end
    end

    return DataFrame("m" => sizes[:], "time" => times[:])
end


function fptas(epsilon, M, n_markets)
    sizes = zeros(Int64, n_markets, length(M))
    times = zeros(Float64, n_markets, length(M))

    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for (j, m) in enumerate(M)
            f, t, g, H = make_correlated_market(m)
            sizes[i, j] = m
            times[i, j] = 1000 * minimum(@elapsed optimalportfolio_fptas(f, t, g, H, epsilon) for r in 1:n_reps)
        end
    end

    return DataFrame("m" => sizes[:], "time" => times[:])
end



function table1()
    printheader("Benchmark 1: Homogeneous-cost algorithms")
    M = fullscale ? [5, 50, 500, 5000] : [5, 10]
    n_markets = fullscale ? 50 : 10

    println("Timing with dictionary")
    grouped_het_dict = groupby(samecosts(:dict, M, n_markets), :m)

    println("Timing with heap")
    grouped_het_heap = groupby(samecosts(:heap, M, n_markets), :m)

    return DataFrame("m" => M,
        "dict_time_mean" => combine(grouped_het_dict, :time => mean)[!, :time_mean],
        "dict_time_std" => combine(grouped_het_dict, :time => std)[!, :time_std],
        "heap_time_mean" => combine(grouped_het_heap, :time => mean)[!, :time_mean],
        "heap_time_std" => combine(grouped_het_heap, :time => std)[!, :time_std])
end


function table2()
    printheader("Benchmark 2: Heterogeneous-cost algorithms")
    M = fullscale ? [5, 50, 500] : [5, 10]
    n_markets = fullscale ? 50 : 10
    epsilons = [0.5, 0.05]

    println("Timing dynamic program")
    dp_times = groupby(dp_H(M, n_markets), :m)
    fptas_times = GroupedDataFrame[]
    for e in epsilons
        println("Timing FPTAS with ε = $e")
        push!(fptas_times, groupby(fptas(e, M, n_markets), :m))
    end

    return DataFrame("m" => M,
        "DP" => combine(dp_times, :time => mean)[!, :time_mean],
        ["fptas_$(epsilons[i])_$f" => combine(fptas_times[i], :time => f)[!, Symbol("time_$f")] for f in [mean, std], i in 1:length(epsilons)]...)
end

println()
display(table1())
println("\n")
display(table2())