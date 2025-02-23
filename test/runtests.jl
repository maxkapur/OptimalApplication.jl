using OptimalApplication
using Test: @test, @testset, @test_throws
using Random: randperm

import Aqua

# Number of markets to generate for each test
const n_markets = 3

@testset verbose = true "OptimalApplication.jl" begin
    @testset verbose = true "Aqua.jl" begin
        Aqua.test_all(OptimalApplication)
    end

    @testset verbose = true "Market I/O" begin
        mkt = SameCostsMarket([0.39, 0.61, 0.03, 0.73, 0.1], [1, 6, 3, 7, 2], 3)
        @test valuation(Int[], mkt) == 0
        @test valuation([1, 3, 5], mkt) ≈ 0.62447

        mkt = SameCostsMarket([0.39, 0.61, 0.03, 0.73, 0.1], Float64[1, 6, 3, 7, 2], 3)
        @test valuation(Int[], mkt) == 0
        @test valuation([1, 3, 5], mkt) ≈ 0.62447

        mkt = VariedCostsMarket(
            [0.39, 0.61, 0.03, 0.73, 0.1],
            [1, 6, 3, 7, 2],
            ones(Int, 5),
            3,
        )
        @test valuation(Int[], mkt) == 0
        @test valuation([1, 3, 5], mkt) ≈ 0.62447
    end

    @testset verbose = true "Same app. costs" begin
        @testset verbose = true "t int" begin
            m = 20

            for _ = 1:n_markets
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

        @testset verbose = true "t float" begin
            m = 20

            for _ = 1:n_markets
                mkt = SameCostsMarket(rand(16), 12 * rand(16), m ÷ 4)

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
    end

    @testset verbose = true "Varied app. costs" begin
        m = 16

        @testset "Exact algorithms" begin
            for _ = 1:n_markets
                mkt = VariedCostsMarket(m)

                X_enum, vX_enum = optimalportfolio_enumerate(mkt)
                sort!(X_enum)
                X_dp, vX_dp = optimalportfolio_dynamicprogram(mkt)
                X_bnb, vX_bnb = optimalportfolio_branchbound(mkt)

                @test X_enum == sort(X_dp)
                @test vX_enum ≈ vX_dp
                @test X_enum == sort(X_bnb)
                @test vX_enum ≈ vX_bnb
            end
        end

        δ = 0.05

        @testset "Slow DP" begin
            for _ = 1:n_markets
                mkt = VariedCostsMarket(m)

                W, vW = optimalportfolio_dynamicprogram_slow(mkt, δ)
                Y, vY = optimalportfolio_dynamicprogram(mkt)

                @test vW ≥ vY - δ
                @test vY == valuation(Y, mkt)
            end
        end

        m = 50
        ε = 0.1

        @testset "FPTAS" begin
            for _ = 1:n_markets
                mkt = VariedCostsMarket(m)

                W, vW = optimalportfolio_fptas(mkt, ε)
                Y, vY = optimalportfolio_dynamicprogram(mkt)

                @test vW / vY ≥ 1 - ε
                @test vY == valuation(Y, mkt)
            end
        end

    end

    # Verify that permuting the problem data doesn't change the optimum
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
            @test v ≈ valuation(optimalportfolio_enumerate(mkt1)[1], mkt1)
            @test v ≈ valuation(optimalportfolio_enumerate(mkt2)[1], mkt2)

            @test v ≤ 2 * valuation(optimalportfolio_fptas(mkt1, 0.5)[1], mkt1)
            @test v ≤ 2 * valuation(optimalportfolio_fptas(mkt2, 0.5)[1], mkt2)
            @test v ≤ 2 * valuation(optimalportfolio_fptas(mkt3, 0.5)[1], mkt3)

            @test v ≤
                  valuation(optimalportfolio_dynamicprogram_slow(mkt1, 0.1)[1], mkt1) + 0.1
            @test v ≤
                  valuation(optimalportfolio_dynamicprogram_slow(mkt2, 0.1)[1], mkt2) + 0.1
            @test v ≤
                  valuation(optimalportfolio_dynamicprogram_slow(mkt3, 0.1)[1], mkt3) + 0.1
        end
    end

    @testset verbose = true "Edge cases" begin
        @testset verbose = true "Same costs" begin
            # h > m should have all schools
            mkt = SameCostsMarket(rand(4), rand(4), 6)

            @test 1:4 == sort(applicationorder_heap(mkt)[1])
            @test 1:4 == sort(applicationorder_list(mkt)[1])
            @test 1:4 == sort(optimalportfolio_enumerate(mkt)[1])
        end

        @testset verbose = true "Het. costs" begin
            # H ≥ sum(g) should have all schools
            mkt = VariedCostsMarket(fill(0.5, 4), fill(2, 4), fill(4, 4), 16)

            @test 1:4 == sort(optimalportfolio_branchbound(mkt)[1])
            @test 1:4 == sort(optimalportfolio_dynamicprogram(mkt)[1])
            @test 1:4 == sort(optimalportfolio_dynamicprogram_slow(mkt, 0.1)[1])
            @test 1:4 == sort(optimalportfolio_fptas(mkt, 0.05)[1])
            @test 1:4 == sort(optimalportfolio_enumerate(mkt)[1])
            @test 1:4 == sort(optimalportfolio_greedy(mkt)[1])
            @test 1:4 == sort(optimalportfolio_simulatedannealing(mkt; nit = 4)[1])

            # g[j] > H should be excluded
            mkt = VariedCostsMarket([0.5, 0.5, 1.0], [1, 2, 100], [3, 3, 5], 4)

            @test [2] == optimalportfolio_branchbound(mkt)[1]
            @test [2] == optimalportfolio_dynamicprogram(mkt)[1]
            @test [2] == optimalportfolio_dynamicprogram_slow(mkt, 0.1)[1]
            @test [2] == optimalportfolio_fptas(mkt, 0.25)[1]
            @test [2] == optimalportfolio_enumerate(mkt)[1]
            @test [2] == sort(optimalportfolio_greedy(mkt)[1])
            @test [2] == sort(optimalportfolio_simulatedannealing(mkt; nit = 4)[1])
        end
    end

    @testset verbose = true "Heuristic algorithms" begin
        mkt = VariedCostsMarket([0.34, 0.8, 0.1, 0.5], [9, 2, 2, 4], [3, 4, 4, 3], 7)
        X, v = optimalportfolio_greedy(mkt)
        @test v == valuation(X, mkt)
        @test sort(X) == [1, 4]
        Y, w = optimalportfolio_simulatedannealing(mkt; nit = 10)
        @test w == valuation(Y, mkt)
        @test w ≥ v

        mkt =
            VariedCostsMarket([0.131, 0.136, 0.09, 0.016], [2, 3, 9, 54], [9, 5, 5, 9], 14)
        X, v = optimalportfolio_greedy(mkt)
        @test v == valuation(X, mkt)
        @test sort(X) == [3, 4]
        Y, w = optimalportfolio_simulatedannealing(mkt; nit = 10)
        @test w == valuation(Y, mkt)
        @test w ≥ v

        m = 20
        for _ = 1:n_markets
            mkt = VariedCostsMarket(m)

            X, v = optimalportfolio_greedy(mkt)
            @test v == valuation(X, mkt)
            Y, w = optimalportfolio_simulatedannealing(mkt; nit = 50)
            @test w == valuation(Y, mkt)
            @test w ≥ v
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
        @test !isempty(optimalportfolio_greedy(mkt)[1])
        @test !isempty(optimalportfolio_simulatedannealing(mkt; nit = 10)[1])
    end

    @testset verbose = true "Bad markets" begin
        # Dim mismatch
        f = [0.1, 0.1]
        t = [1, 2, 4]
        g = [2, 2]
        H = 3

        @test_throws DimensionMismatch Market(f, t, 1)
        @test_throws DimensionMismatch Market(f, t, g, H)

        # f not in (0, 1]
        f = [5.0, 1.0]
        t = [4, 7]

        @test_throws DomainError Market(f, t, 1)
        @test_throws DomainError Market(f, t, g, H)

        # t negative
        f = [5.0, 1.0]
        t = [4, -7]

        @test_throws DomainError Market(f, t, 1)
        @test_throws DomainError Market(f, t, g, H)
    end
end
