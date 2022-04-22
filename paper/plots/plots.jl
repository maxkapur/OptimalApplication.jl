using OptimalApplication
using Random
using CairoMakie
using DataFrames
using Colors

# const accentcolor = (107,142,35) ./ 255
const accentcolor = (210, 10, 60) ./ 255

#   Comparative statics example
#   =============================

function comparativestatics()
    mkt = Market([0.39, 0.33, 0.24, 0.24, 0.05, 0.03, 0.1, 0.12],
        collect(range(start=200, step=50, length=8)),
        8)

    apporder, v = applicationorder_list(mkt, verbose=true)

    invpermapporder = invperm(apporder)

    println(apporder)

    df = DataFrame("f" => mkt.f, "t" => mkt.t, "order" => invpermapporder, "v" => v[invpermapporder])
    #=
     Row │ f        t        order  v       
         │ Float64  Float64  Int64  Float64 
    ─────┼──────────────────────────────────
       1 │    0.39    200      4  230.047
       2 │    0.33    250      2  146.7
       3 │    0.24    300      6  281.513
       4 │    0.24    350      1   84.0
       5 │    0.05    400      7  288.778
       6 │    0.03    450      8  294.106
       7 │    0.1     500      5  257.643
       8 │    0.12    550      3  195.096
    =#
    display(df)

    fig1 = Figure(resolution=(600, 450))

    ax = Axis(fig1[1, 1], xticks=0:8, xlabel="h", ylabel="v*")

    scatterlines!(ax,
        0:8,
        vcat(0, v),
        color=RGB(accentcolor...),
        marker=:diamond,
        markersize=12,
        linewidth=2,
    )

    return fig1
end

#   Sample random market
#   ======================


function samplerandom()
    m = 64
    mkt = VariedCostsMarket(m)

    X, v = optimalportfolio_dynamicprogram(mkt)
    X[:] = invperm(mkt.perm)[X]

    nX = setdiff(1:m, X)

    # GR plots
    # scal = 1.5
    # pl = plot(size = (500, 500), xlabel = "fⱼ", ylabel = "tⱼ", legend = :topright)
    # scatter!(pl, mkt.f[nX], mkt.t[nX], ms = scal * sqrt.(mkt.g[nX]), c = :gray25, msc = :auto, ma = 0.5, label = nothing)
    # scatter!(pl, mkt.f[X], mkt.t[X], m = :utriangle, ms = scal * sqrt.(mkt.g[X]), c = accentcolor, msc = :auto, ma = 0.85, label = "Apply when H = $(mkt.H)")
    # annotate!(pl, [(maximum(mkt.f), 0.9 * maximum(mkt.t), text("(Marker area: gⱼ)  ", 9, :right))])

    scal = 3
    fig2 = Figure(resolution=(600, 600))

    ax = Axis(fig2[1, 1], xlabel="fⱼ", ylabel="tⱼ", xlabelsize=18, ylabelsize=18)

    scatter!(mkt.f[nX], mkt.t[nX], markersize=scal * sqrt.(mkt.g[nX]), color=RGBA(0.25, 0.25, 0.25, 0.8)) #:gray25, strokealpha=0.5)
    scatter!(mkt.f[X], mkt.t[X], marker=:utriangle, markersize=1.6 * scal * sqrt.(mkt.g[X]), color=RGBA(accentcolor..., 0.7), label="Apply when H = $(mkt.H)")

    axislegend(ax, position=:rt, labelsize=16)
    text!("(Marker area: gⱼ)", position=(maximum(mkt.f), 0.9 * maximum(mkt.t)), align=(:right, :baseline), textsize=16)

    return fig2
end

#   1D simulated annealing accuracy
#   =================================

function accuracy_simulatedannealing()
    n_markets = 5000
    mkt_size = 64

    mkts = VariedCostsMarket[VariedCostsMarket(mkt_size) for _ in 1:n_markets]

    v_siman = zeros(n_markets)
    v_exact = zeros(n_markets)

    Threads.@threads for i in 1:n_markets
        v_siman[i] = optimalportfolio_simulatedannealing(mkts[i])[2]
        v_exact[i] = optimalportfolio_dynamicprogram(mkts[i])[2]
    end

    rats = v_siman ./ v_exact

    fig3 = Figure(
        resolution=(600, 450),
        figure_padding=(15, 25, 15, 15)
    )

    ax = Axis(
        fig3[1, 1],
        xlabel="accuracy ratio",
        ylabel="count",
        xticks=0:0.01:1.0,
    )

    hist!(ax,
        rats,
        bins=range(minimum(rats), 1, length=20),
        color=RGB(accentcolor...),
        # bar_labels=:values,
        # label_formatter=x -> round(Int, x),
        # label_size=12,
        # strokewidth=0.5,
        # strokecolor=RGB(accentcolor...),
    )

    xlims!(ax, minimum(rats), 1)
    # xticks!(fig3, xtickrange=(0.9, 1), xticklabels=collect(0.9:0.01:1))
    ylims!(ax, 0, nothing)

    return fig3
end


#   2D simulated annealing accuracy
#   =================================

function accuracy_simulatedannealing_2d()
    n_markets = 500

    sizerangelog2 = (3, 11)

    mkt_sizes = round.(
        Int,
        2.0 .^ (
            sizerangelog2[1] .+ (sizerangelog2[2] - sizerangelog2[1]) * rand(n_markets)
        )
    )

    mkts = VariedCostsMarket[VariedCostsMarket(s) for s in mkt_sizes]

    v_siman = zeros(n_markets)
    v_exact = zeros(n_markets)

    Threads.@threads for i in 1:n_markets
        v_siman[i] = optimalportfolio_simulatedannealing(mkts[i])[2]
        v_exact[i] = optimalportfolio_dynamicprogram(mkts[i])[2]
    end

    rats = v_siman ./ v_exact

    fig4 = Figure(
        resolution=(600, 450),
        figure_padding=(15, 25, 15, 15)
    )

    ax = Axis(
        fig4[1, 1],
        xlabel="m",
        ylabel="accuracy ratio",
        xscale=log2,
        xticks=2 .^ range(sizerangelog2...),
    )

    scatter!(ax,
        mkt_sizes,
        rats,
        bins=range(minimum(rats), 1, length=20),
        color=RGBA(accentcolor..., 0.7),
        # bar_labels=:values,
        # label_formatter=x -> round(Int, x),
        # label_size=12,
        markersize=5,
    )

    # ylims!(ax, minimum(rats), 1)
    # # xticks!(fig3, xtickrange=(0.9, 1), xticklabels=collect(0.9:0.01:1))
    # ylims!(ax, 0, nothing)

    return fig4
end


# fig1 = comparativestatics()
# save("paper/plots/h_v-example.pdf", fig1)
# save("paper/plots/h_v-example.png", fig1)

# fig2 = samplerandom()
# save("paper/plots/samplemarket.pdf", fig2)
# save("paper/plots/samplemarket.png", fig2)

# fig3 = accuracy_simulatedannealing()
# save("paper/plots/accuracy_simulatedannealing.pdf", fig3)
# save("paper/plots/accuracy_simulatedannealing.png", fig3)

# fig4 = accuracy_simulatedannealing_2d()
# save("paper/plots/accuracy_simulatedannealing.pdf", fig4)
# save("paper/plots/accuracy_simulatedannealing.png", fig4)
