using OptimalApplication
using DataFrames
using PrettyTables
using Plots
using Statistics: median, quantile

function plot_params()
    Ts = 2.0 .^ (-4:2:4)
    rs = (2.0 .^ (-5:0))

    n_markets = 1000

    sizerangelog2 = (3, 9)

    sizes =
        round.(
            Int,
            2.0 .^
            (sizerangelog2[1] .+ (sizerangelog2[2] - sizerangelog2[1]) * rand(n_markets)),
        )

    rats = zeros(length(rs), length(Ts), n_markets)

    @time for k = 1:n_markets
        mkt = VariedCostsMarket(sizes[k])
        opt = optimalportfolio_dynamicprogram(mkt)[2]

        Threads.@threads for j = 1:length(Ts)
            for i = 1:length(rs)
                rats[i, j, k] =
                    optimalportfolio_simulatedannealing(mkt; temp = Ts[j], red = rs[i])[2]
            end
        end

        rats[:, :, k] /= optimalportfolio_dynamicprogram(mkt)[2]
    end

    xts = 2 .^ range(sizerangelog2...)

    pl = plot(
        layout = grid(length(rs), length(Ts)),
        link = :all,
        xscale = :log10,
        xticks = (xts, xts),
        legend = nothing,
        size = (1200, 1200),
        xlim = 2 .^ (sizerangelog2 .+ (-0.5, 0.5)),
    )


    avgs = zeros(length(rs), length(Ts))
    for j = 1:length(Ts), i = 1:length(rs)
        avgs[i, j] = quantile(rats[i, j, :], 0.1)
    end
    # [sum(rats[i, j, :]) / length(rats[i, j, :]) for i in 1:length(rs), j in 1:length(Ts)]
    best_avg = maximum(avgs)
    @show best_avg

    for j = 1:length(Ts), i = 1:length(rs)
        hline!(pl[i, j], [avgs[i, j]], c = :navy)
        scatter!(
            pl[i, j],
            sizes,
            rats[i, j, :],
            c = avgs[i, j] ≥ best_avg ? :olivedrab : :crimson,
            msc = :auto,
            alpha = 0.7,
            title = "T = $(round(Ts[j], digits=2)), r = $(round(rs[i], digits=4))",
            titlefontsize = 12,
        )
        if i < length(rs)
            plot!(pl[i, j], xticks = (xts, fill("", length(xts))))
        end
        if j > 1
            plot!(
                pl[i, j],
                yticks = (Plots.get_ticks(pl[1, 1], :y)[1], fill("", length(xts))),
            )
        end
    end

    return pl, rats, sizes
end

pl, rats, sizes = plot_params()
display(pl)
# savefig("benchmark/annealingparams-plot.png")
# savefig("benchmark/annealingparams-plot.pdf")
