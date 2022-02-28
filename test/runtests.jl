using OptimalApplication
using Test

# Number of markets to generate for each test
n_markets = 3

function make_correlated_market(m)
    φ = 0.7
    t = rand(m) |> sort
    f = 1 .- (φ * t + (1 - φ) * rand(m))
    return f, t
end

@testset verbose = true "Homogeneous application costs" begin
    m = 20
    h = m ÷ 2

    for _ in 1:n_markets
        f, t = make_correlated_market(m)

        X, vX = optimalportfolio_enumerate(f, t, h)
        sort!(X)
        W, vW = applicationorder(f, t, h; datastructure = :heap)
        Y, vY = applicationorder(f, t, h; datastructure = :dict)
        @test X == sort(W)
        @test vX ≈ vW[h]
        @test X == sort(Y)
        @test vX ≈ vY[h]
    end
end

@testset verbose = true "Heterogeneous application costs" begin
    m = 20

    @testset "Exact algorithms" begin
        for _ in 1:n_markets
            f, t = make_correlated_market(m)
            g = rand(1:5, m)
            H = sum(g) ÷ 2

            X, vX = optimalportfolio_enumerate(f, t, g, H)
            sort!(X)
            W, VW = optimalportfolio_valuationtable(f, t, g, H)
            Y, vY = optimalportfolio_dynamicprogram(f, t, g, H)
            # Without memoization
            Z, vZ = optimalportfolio_dynamicprogram(f, t, g, H, false)

            @test X == sort(W)
            @test vX ≈ VW[m, H]
            @test X == sort(Y)
            @test vX ≈ vY
            @test X == sort(Z)
            @test vX ≈ vZ
        end
    end

    m = 100
    ε = 0.1

    @testset "FPTAS" begin
        for _ in 1:n_markets
            f, t = make_correlated_market(m)
            t = round.(Int, 100 * t)
            g = rand(1:5, m)
            H = sum(g) ÷ 2

            W, vW = optimalportfolio_fptas(f, t, g, H, ε)
            Y, vY = optimalportfolio_dynamicprogram(f, t, g, H)

            @test vW / vY ≥ 1 - ε
        end
    end
end
