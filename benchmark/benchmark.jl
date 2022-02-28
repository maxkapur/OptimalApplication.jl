using OptimalApplication
using DataFrames
using Random
using StatsBase
using Base.Threads

# An all-day benchmark; tweak parameters of run with caution.

# Number of times to repeat each computation, where min of these is reported as time
n_reps = 3

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
    sizes = Int64[]
    times = Float64[]

    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for m in M
            f, t, = make_correlated_market(m)
            push!(sizes, m)
            push!(times, minimum(@elapsed applicationorder(f, t; datastructure = ds) for r in 1:n_reps))
        end
    end

    return DataFrame("m" => sizes, "time" => times)
end


function dp_H(M, n_markets)
    sizes = Int64[]
    times = Float64[]

    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for m in M
            f, t, g, H = make_correlated_market(m)
            push!(sizes, m)
            push!(times, minimum(@elapsed optimalportfolio_dynamicprogram(f, t, g, H) for r in 1:n_reps))
        end
    end

    return DataFrame("m" => sizes, "time" => times)
end


function fptas(epsilon, M, n_markets)
    sizes = Int64[]
    times = Float64[]

    @threads for i in 1:n_markets
        println("  i = $i of $n_markets")
        for m in M
            f, t, g, H = make_correlated_market(m)
            push!(sizes, m)
            push!(times, minimum(@elapsed optimalportfolio_fptas(f, t, g, H, epsilon) for r in 1:n_reps))
        end
    end

    return DataFrame("m" => sizes, "time" => times)
end



function table1()
    println("== Homogeneous cost algorithms ==")
    M = [5, 50, 500, 5000]
    n_markets = 50

    println("Timing with dictionary")
    het_dict = samecosts(:dict, M, n_markets)

    println("Timing with heap")
    het_heap = samecosts(:heap, M, n_markets)

    return DataFrame("m" => M,
        "dict_time" => combine(groupby(het_dict, :m), :time => mean)[!, :time_mean],
        "heap_time" => combine(groupby(het_heap, :m), :time => mean)[!, :time_mean])
end


function table2()
    println("== Heterogeneous cost algorithms ==")
    M = [10, 100, 1000]
    n_markets = 50
    epsilons = [0.5, 0.05]

    println("Timing dynamic program")
    dp_times = dp_H(M, n_markets)
    fptas_times = DataFrame[]
    for e in epsilons
        println("Timing FPTAS with ε = $e")
        push!(fptas_times, fptas(e, M, n_markets))
    end

    @show n_markets

    return DataFrame("m" => M,
        "DP" => combine(groupby(dp_times, :m), :time => mean)[!, :time_mean],
        ["fptas_$e" => combine(groupby(fptas_times[i], :m), :time => mean)[!, :time_mean] for (i, e) in enumerate(epsilons)]...)
end

display(table1())
display(table2())
