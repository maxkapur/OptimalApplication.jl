using OptimalApplication
using Test
using Random

# Number of markets to generate for each test
n_markets = 3

function make_correlated_market(m)
    A = 10
    t = ceil.(Int, 10 * randexp(m))
    sort!(t)
    f = inv.(t .+ 10 * rand(m))
    g = rand(5:10, m)
    H = sum(g) ÷ 2
    return f, t, g, H
end

@testset verbose = true "OptimalApplication.jl" begin
    @testset verbose = true "Same app. costs" begin
        m = 20
        h = m ÷ 2

        for _ in 1:n_markets
            f, t, = make_correlated_market(m)
            if rand() > 0.5
                t = Float64.(t)
            end

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

    @testset verbose = true "Varied app. costs" begin
        m = 20

        @testset "Exact algorithms" begin
            for _ in 1:n_markets
                f, t, g, H = make_correlated_market(m)
                if rand() > 0.5
                    t = Float64.(t)
                end
            
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
                f, t, g, H = make_correlated_market(m)

                W, vW = optimalportfolio_fptas(f, t, g, H, ε)
                Y, vY = optimalportfolio_dynamicprogram(f, t, g, H)

                @test vW / vY ≥ 1 - ε
            end
        end
    end
end
