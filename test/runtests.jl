using OptimalApplication
using Test
using Random

# "t not sorted" warnings are OK. 

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

    # These tests verify that permuting the problem data doesn't
    # change the optimum
    @testset verbose = true "Permutation" begin
        @testset verbose = true "Same costs" begin
            m = 6
            f = rand(m)
            t = rand(m)
            h = m ÷ 2

            perm = randperm(m)
            sp = sortperm(t)

            # Randomly sorted market
            mkt1 = SameCostsMarket(f, t, h)
            # A random permutation thereof
            mkt2 = SameCostsMarket(f[perm], t[perm], h)
            # Now sorted by t
            mkt3 = SameCostsMarket(f[sp], t[sp], h)

            # A random portfolio
            x = rand(Bool, m)
            @test valuation((1:m)[x], mkt1) ≈ valuation((1:m)[x[perm]], mkt2)
            @test valuation((1:m)[x], mkt1) ≈ valuation((1:m)[x[sp]], mkt3)

            # Optimal valuations all agree
            x, v = optimalportfolio_enumerate(mkt3)
            @test v ≈ valuation(applicationorder_list(mkt1)[1], mkt1)
            @test v ≈ valuation(applicationorder_list(mkt2)[1], mkt2)
            @test v ≈ valuation(applicationorder_list(mkt3)[1], mkt3)
            @test v ≈ valuation(applicationorder_heap(mkt1)[1], mkt1)
            @test v ≈ valuation(applicationorder_heap(mkt2)[1], mkt2)
            @test v ≈ valuation(applicationorder_heap(mkt3)[1], mkt3)
            @test v ≈ valuation(optimalportfolio_enumerate(mkt1)[1], mkt1)
            @test v ≈ valuation(optimalportfolio_enumerate(mkt2)[1], mkt2)
        end

        @testset verbose = true "Het. costs" begin
            m = 6
            f = rand(m)
            t = rand(50:60, m)
            g = rand(5:10, m)
            H = sum(g) ÷ 2

            perm = randperm(m)
            sp = sortperm(t)

            # Randomly sorted market
            mkt1 = VariedCostsMarket(f, t, g, H)
            # A random permutation thereof
            mkt2 = VariedCostsMarket(f[perm], t[perm], g[perm], H)
            # Now sorted by t
            mkt3 = VariedCostsMarket(f[sp], t[sp], g[sp], H)

            # A random portfolio
            x = rand(Bool, m)
            @test valuation((1:m)[x], mkt1) ≈ valuation((1:m)[x[perm]], mkt2)
            @test valuation((1:m)[x], mkt1) ≈ valuation((1:m)[x[sp]], mkt3)

            # Optimal valuations all agree
            x, v = optimalportfolio_enumerate(mkt3)
            @test v ≈ valuation(optimalportfolio_branchbound(mkt1)[1], mkt1)
            @test v ≈ valuation(optimalportfolio_branchbound(mkt2)[1], mkt2)
            @test v ≈ valuation(optimalportfolio_branchbound(mkt3)[1], mkt3)
            @test v ≈ valuation(optimalportfolio_dynamicprogram(mkt1)[1], mkt1)
            @test v ≈ valuation(optimalportfolio_dynamicprogram(mkt2)[1], mkt2)
            @test v ≈ valuation(optimalportfolio_dynamicprogram(mkt3)[1], mkt3)
            @test v ≈ valuation(optimalportfolio_valuationtable(mkt1)[1], mkt1)
            @test v ≈ valuation(optimalportfolio_valuationtable(mkt2)[1], mkt2)
            @test v ≈ valuation(optimalportfolio_valuationtable(mkt3)[1], mkt3)
            @test v ≈ valuation(optimalportfolio_enumerate(mkt1)[1], mkt1)
            @test v ≈ valuation(optimalportfolio_enumerate(mkt2)[1], mkt2)

            @test v ≤ 2 * valuation(optimalportfolio_fptas(mkt1, 0.5)[1], mkt1)
            @test v ≤ 2 * valuation(optimalportfolio_fptas(mkt2, 0.5)[1], mkt2)
            @test v ≤ 2 * valuation(optimalportfolio_fptas(mkt3, 0.5)[1], mkt3)
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
