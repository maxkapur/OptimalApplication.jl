function iscoherentmarket(f::Vector{Float64}, t::Vector{Int})
    @assert length(f) == length(t)
    @assert all(0 .< f .≤ 1)
    @assert issorted(t)
end


function iscoherentmarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int})
    @assert length(f) == length(t)
    @assert length(f) == length(g)
    @assert all(0 .< f .≤ 1)
    @assert issorted(t)
end


function isnontrivialmarket(f::Vector{Float64}, t::Vector{Int}, h::Int)
    iscoherentmarket(f, t)
    @assert 0 < h ≤ length(t)
end


function isnontrivialmarket(
    f::Vector{Float64},
    t::Vector{Int},
    g::Vector{Int},
    H::Int)
    iscoherentmarket(f, t, g)
    @assert 0 < H ≤ sum(g)
    @assert all(0 .< g .≤ H)
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
struct SameCostsMarket
    m::Int
    f::Vector{Float64}
    t::Vector{Int}
    h::Int
    ft::Vector{Float64}      # = f .* t
    omf::Vector{Float64}     # = 1 .- f
    perm::Vector{Int}

    """
        SameCostsMarket(f, t, h)

    Perform preliminary helper calculations for `SameCostsMarket` defined by
    admissions probabilities `f`, utility values `t`, and application limit `h`
    and return the market object.
    """
    function SameCostsMarket(f::Vector{Float64}, t::Vector{Int}, h::Int)
        isnontrivialmarket(f, t, h)
        m = length(f)

        # We are current asserting that t is sorted; a placeholder for later work
        perm = sortperm(t)

        return new(m, f, t, h, f .* t, 1 .- f, perm)
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
struct VariedCostsMarket
    m::Int
    f::Vector{Float64}
    t::Vector{Int}
    g::Vector{Int}
    H::Real
    ft::Vector{Float64}      # = f .* t
    omf::Vector{Float64}     # = 1 .- f
    perm::Vector{Int}

    """
        VariedCostsMarket(f, t, g, H)

    Perform preliminary helper calculations for `SameCostsMarket` defined by
    admissions probabilities `f`, utility values `t`, application costs `g`,
    and application budget `H` and return the market object.
    """
    function VariedCostsMarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int}, H::Int)
        isnontrivialmarket(f, t, g, H)
        m = length(f)
        perm = sortperm(t)

        return new(m, f, t, g, H, f .* t, 1 .- f, perm)
    end
end

"""
    Market(f, t, h)

Perform preliminary helper calculations for `SameCostsMarket` defined by
admissions probabilities `f`, utility values `t`, and application limit `h`
and return the market object.
"""
function Market(f::Vector{Float64}, t::Vector{Int}, h::Int)
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
function valuation(X::Vector{Int}, mkt::Union{SameCostsMarket,VariedCostsMarket})
    isempty(X) && return 0.0

    sort!(X)
    h = length(X)

    if h > 1
        res = 0.0
        cp = reverse(cumprod(reverse(mkt.omf[X[2:end]])))

        for j in 1:h-1
            res += mkt.ft[X[j]] * cp[j]
        end

        res += mkt.ft[X[end]]

        return res
    else
        return mkt.ft[X[1]]
    end
end
