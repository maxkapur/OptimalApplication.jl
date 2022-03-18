using OptimalApplication
using Test
using Random

# Number of markets to generate for each test
const n_markets = 3

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
    f, t, g, H = make_correlated_market(m)
    return SameCostsMarket(f, t, m ÷ 2)
end

@testset verbose = true "OptimalApplication.jl" begin
    @testset verbose = true "Same app. costs" begin
        m = 20

        for _ in 1:n_markets
            mkt = randSCM(m)

            X, vX = optimalportfolio_enumerate(mkt)
            sort!(X)
            W, vW = applicationorder_heap(mkt)
            Y, vY = applicationorder_list(mkt)
            @test X == sort(W)
            @test vX ≈ last(vW)
            @test X == sort(Y)
            @test vX ≈ last(vY)
        end
    end

    @testset verbose = true "Varied app. costs" begin
        m = 10

        @testset "Exact algorithms" begin
            for _ in 1:n_markets
                mkt = randVCM(m)

                X, vX = optimalportfolio_enumerate(mkt)
                sort!(X)
                W, VW = optimalportfolio_valuationtable(mkt)
                Y, vY = optimalportfolio_dynamicprogram(mkt)
                B, vB = optimalportfolio_branchbound(mkt)

                @test X == sort(W)
                @test vX ≈ VW[mkt.m, mkt.H]
                @test X == sort(Y)
                @test vX ≈ vY
                @test X == sort(B)
                @test vX ≈ vB
            end
        end

        m = 50
        ε = 0.1

        @testset "FPTAS" begin
            for _ in 1:n_markets
                mkt = randVCM(m)

                W, vW = optimalportfolio_fptas(mkt, ε)
                Y, vY = optimalportfolio_dynamicprogram(mkt)

                @test vW / vY ≥ 1 - ε
            end
        end
    end

    @testset verbose = true "Large problems" begin
        mkt = randSCM(5000)

        X, V = applicationorder_list(mkt)
        @test !isempty(X)

        mkt = randVCM(500)

        W, vW = optimalportfolio_fptas(mkt, 0.5)
        Y, vY = optimalportfolio_dynamicprogram(mkt)
        @test vW / vY ≥ 0.5
    end

    @testset verbose = true "Bad markets" begin
        # Dim mismatch
        f = [0.1, 0.1]
        t = [1, 2, 4]
        g = [2, 2]
        H = 3

        @test_throws AssertionError Market(f, t, 1)
        @test_throws AssertionError Market(f, t, g, H)

        # t not integer
        f = [0.1, 0.1]
        t = [4.0, 7.0]
        g = [2, 2]
        H = 3

        @test_throws MethodError Market(f, t, 1)
        @test_throws MethodError Market(f, t, g, H)

        # t not sorted
        f = [0.1, 0.1]
        t = [7, 4]
        g = [2, 2]
        H = 3

        @test_throws AssertionError Market(f, t, 1)
        @test_throws AssertionError Market(f, t, g, H)

        # f not in (0, 1]
        f = [5.0, 1.0]
        t = [4, 7]

        @test_throws AssertionError Market(f, t, 1)
        @test_throws AssertionError Market(f, t, g, H)

        # Some g[j] > H
        f = [0.1, 0.1]
        t = [4, 7]
        g = [10, 10]
        H = 5

        @test_throws AssertionError Market(f, t, g, H)
    end
end