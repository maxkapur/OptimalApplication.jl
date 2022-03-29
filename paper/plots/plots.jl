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

    nX = setdiff(1:m, X)

    # GR plots
    # scal = 1.5
    # pl = plot(size = (500, 500), xlabel = "fⱼ", ylabel = "tⱼ", legend = :topright)
    # scatter!(pl, mkt.f[nX], mkt.t[nX], ms = scal * sqrt.(mkt.g[nX]), c = :gray25, msc = :auto, ma = 0.5, label = nothing)
    # scatter!(pl, mkt.f[X], mkt.t[X], m = :utriangle, ms = scal * sqrt.(mkt.g[X]), c = accentcolor, msc = :auto, ma = 0.85, label = "Apply when H = $(mkt.H)")
    # annotate!(pl, [(maximum(mkt.f), 0.9 * maximum(mkt.t), text("(Marker area: gⱼ)  ", 9, :right))])

    scal = 3
    fig2 = Figure(resolution=(600, 600))

    ax = Axis(fig2[1, 1], xlabel = "fⱼ", ylabel = "tⱼ", xlabelsize=18, ylabelsize=18)

    scatter!(mkt.f[nX], mkt.t[nX], markersize=scal * sqrt.(mkt.g[nX]), color=RGBA(0.25, 0.25, 0.25, 0.8)) #:gray25, strokealpha=0.5)
    scatter!(mkt.f[X], mkt.t[X], marker=:utriangle, markersize=1.6 * scal * sqrt.(mkt.g[X]), color=RGBA(accentcolor..., 0.85), label="Apply when H = $(mkt.H)")

    axislegend(ax, position=:rt, labelsize=16)
    text!("(Marker area: gⱼ)", position=(maximum(mkt.f), 0.9 * maximum(mkt.t)), align=(:right, :baseline), textsize=16)

    return fig2
end

fig1 = comparativestatics()
fig2 = samplerandom()

save("paper/plots/h_v-example.pdf", fig1)
save("paper/plots/h_v-example.png", fig1)

save("paper/plots/samplemarket.pdf", fig2)
save("paper/plots/samplemarket.png", fig2)
