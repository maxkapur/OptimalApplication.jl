using OptimalApplication
using Test
using Random

# Number of markets to generate for each test
const n_markets = 3

@testset verbose = true "OptimalApplication.jl" begin
    @testset verbose = true "Same app. costs, t int" begin
        m = 20

        for _ in 1:n_markets
            mkt = SameCostsMarket(m)

            X_enum, vX_enum = optimalportfolio_enumerate(mkt)
            sort!(X_enum)
            X_heap, vX_heap = applicationorder_heap(mkt)
            X_list, vX_list = applicationorder_list(mkt)
            @test X_enum == sort(X_heap)
            @test vX_enum ≈ last(vX_heap)
            @test X_enum == sort(X_list)
            @test vX_enum ≈ last(vX_list)
        end
    end

    @testset verbose = true "Same app. costs, t float" begin
        m = 20

        for _ in 1:n_markets
            mkt = SameCostsMarket(rand(16), sort(12 * rand(16)), m ÷ 4)

            X_enum, vX_enum = optimalportfolio_enumerate(mkt)
            sort!(X_enum)
            X_heap, vX_heap = applicationorder_heap(mkt)
            X_list, vX_list = applicationorder_list(mkt)
            @test X_enum == sort(X_heap)
            @test vX_enum ≈ last(vX_heap)
            @test X_enum == sort(X_list)
            @test vX_enum ≈ last(vX_list)
        end
    end

    @testset verbose = true "Varied app. costs" begin
        m = 16

        @testset "Exact algorithms" begin
            for _ in 1:n_markets
                mkt = VariedCostsMarket(m)

                X_enum, vX_enum = optimalportfolio_enumerate(mkt)
                sort!(X_enum)
                X_valtable, VX_valtable = optimalportfolio_valuationtable(mkt)
                X_dp, vX_dp = optimalportfolio_dynamicprogram(mkt)
                X_bnb, vX_bnb = optimalportfolio_branchbound(mkt)

                @test X_enum == sort(X_valtable)
                @test vX_enum ≈ VX_valtable[mkt.m, mkt.H]
                @test X_enum == sort(X_dp)
                @test vX_enum ≈ vX_dp
                @test X_enum == sort(X_bnb)
                @test vX_enum ≈ vX_bnb
            end
        end

        m = 50
        ε = 0.1

        @testset "FPTAS" begin
            for _ in 1:n_markets
                mkt = VariedCostsMarket(m)

                W, vW = optimalportfolio_fptas(mkt, ε)
                Y, vY = optimalportfolio_dynamicprogram(mkt)

                @test vW / vY ≥ 1 - ε
            end
        end
    end

    @testset verbose = true "Large problems" begin
        mkt = SameCostsMarket(5000)

        X, V = applicationorder_list(mkt)
        @test !isempty(X)

        mkt = VariedCostsMarket(500)

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
