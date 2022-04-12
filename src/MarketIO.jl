const UIntTypes = DataType[UInt8, UInt16, UInt32, UInt64]


# TO DO: Use Julian errors rather than assertions here.
function iscoherentmarket(f::Vector{Float64}, t::Vector{<:Real})
    @assert length(f) == length(t)
    @assert all(0 .< f .≤ 1)
    @assert all(0 .≤ t)
end


function iscoherentmarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int})
    @assert length(f) == length(t)
    @assert length(f) == length(g)
    @assert all(0 .< f .≤ 1)
    @assert all(0 .≤ t)
    @assert all(0 .≤ g)
end


"""
Contains information about a college application market with identical
application costs. Fields:

    `m`: Number of schools
    `f`: Vector of admissions probabilities
    `t`: Vector of utilitied values (must be sorted)
    `h`: Number of schools student is allowed to apply to
    `ft = f .* t`
    `omf = 1 .- f`
    `perm = sortperm(perm)`: Currently a placeholder
"""
struct SameCostsMarket{T<:Unsigned, U<:Real}
    m::T
    f::Vector{Float64}
    t::Vector{U}
    h::T
    ft::Vector{Float64}      # = f .* t
    omf::Vector{Float64}     # = 1 .- f
    perm::Vector{T}

    """
        SameCostsMarket(f, t, h)

    Perform preliminary helper calculations for `SameCostsMarket` defined by
    admissions probabilities `f`, utility values `t`, and application limit `h`
    and return the market object.
    """
    function SameCostsMarket(f::Vector{Float64}, t::Vector{U}, h::Integer) where U<:Real
        iscoherentmarket(f, t)
        T = UIntTypes[findfirst(T -> length(f) < typemax(T), UIntTypes)]

        # We are current asserting that t is sorted; a placeholder for later work
        perm = T.(sortperm(t))

        return new{T, U}(
            T(length(f)),
            f[perm],
            t[perm],
            T(min(h, length(f))),
            (f .* t)[perm],
            (1 .- f)[perm],
            perm)
    end

    function SameCostsMarket(m::Integer)
        T = UIntTypes[findfirst(T -> m < typemax(T), UIntTypes)]
        t = ceil.(Int, -10 * log.(rand(m)))
        f = inv.(t .+ 10 * rand(m))
        perm = T.(sortperm(t))
        return new{T,Int}(
            T(m),
            f[perm],
            t[perm],
            T(m ÷ 2),
            (f.*t)[perm],
            (1 .- f)[perm],
            perm)
    end
end


"""
    VariedCostsMarket

Contains information about a college application market with varying
application costs. Fields:

    `m`: Number of schools
    `f`: Vector of admissions probabilities
    `t`: Vector of utilitied values (must be sorted)
    `g`: Vector of application costs
    `H`: Budget to spend on applications
    `ft = f .* t`
    `omf = 1 .- f`
    `perm = sortperm(perm)`: Currently a placeholder
"""
struct VariedCostsMarket{T<:Unsigned}
    m::T
    f::Vector{Float64}
    t::Vector{Int}
    g::Vector{Int}
    H::Int
    ft::Vector{Float64}      # = f .* t
    omf::Vector{Float64}     # = 1 .- f
    perm::Vector{T}

    """
        VariedCostsMarket(f, t, g, H)

    Perform preliminary helper calculations for `SameCostsMarket` defined by
    admissions probabilities `f`, utility values `t`, application costs `g`,
    and application budget `H` and return the market object.
    """
    function VariedCostsMarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int}, H::Int)
        iscoherentmarket(f, t, g)
        T = UIntTypes[findfirst(T -> length(f) < typemax(T), UIntTypes)]
        perm = T.(sortperm(t))

        return new{T}(
            T(length(f)),
            f[perm],
            t[perm],
            g[perm],
            min(H, sum(g)),
            (f.*t)[perm],
            (1 .- f)[perm],
            perm)
    end

    function VariedCostsMarket(m::Integer)
        T = UIntTypes[findfirst(T -> m < typemax(T), UIntTypes)]
        t = ceil.(Int, -10 * log.(rand(m)))
        f = inv.(t .+ 10 * rand(m))
        g = rand(5:10, m)
        H = sum(g) ÷ 2
        perm = T.(sortperm(t))

        return new{T}(
            T(m),
            f[perm],
            t[perm],
            g[perm],
            H,
            (f.*t)[perm],
            (1 .- f)[perm],
            perm)
    end
end

"""
    Market(f, t, h)

Perform preliminary helper calculations for `SameCostsMarket` defined by
admissions probabilities `f`, utility values `t`, and application limit `h`
and return the market object.
"""
function Market(f::Vector{Float64}, t::Vector{<:Real}, h::Int)
    return SameCostsMarket(f, t, h)
end


"""
    Market(f, t, g, H)

Perform preliminary helper calculations for `SameCostsMarket` defined by
admissions probabilities `f`, utility values `t`, application costs `g`,
and application budget `H` and return the market object.
"""
function Market(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int}, H::Int)
    return VariedCostsMarket(f, t, g, H)
end


"""
    valuation(X, mkt)

Return the valuation of the portfolio `X` on the market `mkt`, which may be either 
a `SameCostsMarket` or a `VariedCostsMarket`.
"""
function valuation(
        X::AbstractVector{<:Integer},
        mkt::Union{SameCostsMarket{T},VariedCostsMarket{T}};
        invp::Union{Nothing,AbstractVector{T}}=nothing)::Float64 where T<:Unsigned
    isempty(X) && return 0.0

    X = T.(X)

    # sort!(X)
    if isnothing(invp)
        invp = invperm(mkt.perm)
    end
    X[:] = sort(invp[X])

    h = T(length(X))

    if h > 1
        res = 0.0
        cp = reverse(cumprod(reverse(mkt.omf[X[2:end]])))

        for j in T(1):T(h-1)
            res += mkt.ft[X[j]] * cp[j]
        end

        res += mkt.ft[X[end]]

        return res
    else
        return mkt.ft[X[1]]
    end
end
