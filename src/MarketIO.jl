function iscoherentmarket(f::Vector{<:Real}, t::Vector{<:Real})
    @assert length(f) == length(t)
    @assert all(0 .< f .≤ 1)
    @assert issorted(t)
end


function iscoherentmarket(f::Vector{<:Real}, t::Vector{<:Real}, g::Vector{<:Real})
    @assert length(f) == length(t)
    @assert length(f) == length(g)
    @assert all(0 .< f .≤ 1)
    @assert issorted(t)
end


function isnontrivialmarket(f::Vector{<:Real}, t::Vector{<:Real}, h::Int)
    iscoherentmarket(f, t)
    @assert 0 < h ≤ length(t)
end


function isnontrivialmarket(
    f::Vector{<:Real},
    t::Vector{<:Real},
    g::Vector{<:Real},
    H::Real)
    iscoherentmarket(f, t, g)
    @assert 0 < H ≤ sum(g)
    @assert all(0 .< g .≤ H)
end


"""
Contains information about a college application market with identical
application costs.
"""
struct SameCostsMarket
    m::Integer
    f::Vector{<:Real}
    t::Vector{<:Real}
    h::Integer
    ft::Vector{<:Real}      # = f .* t
    omf::Vector{<:Real}     # = 1 .- f
    perm::Vector{<:Integer}

    """
        SameCostsMarket(f, t, h)

    Perform preliminary helper calculations for `SameCostsMarket` defined by
    admissions probabilities `f`, utility values `t`, and application limit `h`
    and return the market object.
    """
    function SameCostsMarket(f::Vector{<:Real}, t::Vector{<:Real}, h::Integer)
        isnontrivialmarket(f, t, h)
        m = length(f)

        # We are current asserting that t is sorted; a placeholder for later work
        perm = sortperm(t)

        return new(m, f, t, h, f .* t, 1 .- f, perm)
    end
end


"""
Contains information about a college application market with varying
application costs.
"""
struct VariedCostsMarket
    m::Integer
    f::Vector{<:Real}
    t::Vector{<:Real}
    g::Vector{<:Real}
    H::Real
    ft::Vector{<:Real}      # = f .* t
    omf::Vector{<:Real}     # = 1 .- f
    perm::Vector{<:Integer}

    """
        VariedCostsMarket(f, t, g, H)

    Perform preliminary helper calculations for `SameCostsMarket` defined by
    admissions probabilities `f`, utility values `t`, application costs `g`,
    and application budget `H` and return the market object.
    """
    function VariedCostsMarket(f::Vector{<:Real}, t::Vector{<:Real}, g::Vector{<:Real}, H::Real)
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
function Market(f::Vector{<:Real}, t::Vector{<:Real}, h::Integer)
    return SameCostsMarket(f, t, h)
end


"""
    Market(f, t, g, H)

Perform preliminary helper calculations for `SameCostsMarket` defined by
admissions probabilities `f`, utility values `t`, application costs `g`,
and application budget `H` and return the market object.
"""
function Market(f::Vector{<:Real}, t::Vector{<:Real}, g::Vector{<:Real}, H::Real)
    return VariedCostsMarket(f, t, g, H)
end


"""
    valuation(X, mkt)

Return the valuation of the portfolio `X` on the market `mkt`, which may be either 
a `SameCostsMarket` or a `VariedCostsMarket`.
"""
function valuation(X::Vector{<:Integer}, mkt::Union{SameCostsMarket,VariedCostsMarket})
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
