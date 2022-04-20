"""
    smallestunsigned(k::Integer)

Identify the smallest unsigned integer type for which `k < typemax(T)`
"""
function smallestunsignedtype(k::Integer)::DataType
    UIntTypes = DataType[UInt8, UInt16, UInt32, UInt64]

    j = findfirst(UIntTypes) do T
        k < typemax(T)
    end

    return UIntTypes[j]
end


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
    SameCostsMarket{T<:Unsigned,U<:Real}

Contains information about a college application market with identical
application costs.

# Constructors

```julia
SameCostsMarket(f::Vector{Float64}, t::Vector{U}, h::Integer) where {U<:Real}
```
Construct the `SameCostsMarket` defined by admissions probabilities `f`,
utility values `t`, and application limit `h`.

```jldoctest
julia> SameCostsMarket([0.3, 0.2, 0.05], [2, 3, 4], 2)
SameCostsMarket{UInt8, Int64}(0x03, [0.3, 0.2, 0.05], [2, 3, 4], 0x02, [0.6, 0.6000000000000001, 0.2], [0.7, 0.8, 0.95], UInt8[0x01, 0x02, 0x03])
```

```julia
SameCostsMarket(m)
```
Generate a random `SameCostsMarket` with `m` schools.


# Internal API

This type contains the following fields:

 - `m`: Number of schools
 - `f`: Vector of admissions probabilities
 - `t`: Vector of utility values (must be sorted)
 - `h`: Number of schools student is allowed to apply to
 - `ft = f .* t`
 - `omf = 1 .- f`
 - `perm`: How the input data were permuted to sort by `t`

`T` is the type of `h` and `m` and eltype of `perm`. Typically `T` is the
smallest `UInt` type for which `m < typemax(T)`. `U` is the eltype of `t`.

Internally, the schools must be indexed so that `t` is sorted ascending. `perm`
allows the original order to be recovered:

```jldoctest
julia> f = rand(3); t = rand(3); mkt = SameCostsMarket(f, t, 2);

julia> t[mkt.perm] == mkt.t
true

julia> t == mkt.t[invperm(mkt.perm)]
true
```
"""
struct SameCostsMarket{T<:Unsigned,U<:Real}
    m::T
    f::Vector{Float64}
    t::Vector{U}
    h::T
    ft::Vector{Float64}      # = f .* t
    omf::Vector{Float64}     # = 1 .- f
    perm::Vector{T}

    function SameCostsMarket(f::Vector{Float64}, t::Vector{U}, h::Integer) where {U<:Real}
        iscoherentmarket(f, t)
        T = smallestunsignedtype(length(f))
        perm = T.(sortperm(t))

        return new{T,U}(
            T(length(f)),
            f[perm],
            t[perm],
            T(min(h, length(f))),
            (f.*t)[perm],
            (1 .- f)[perm],
            perm)
    end

    function SameCostsMarket(m::Integer)
        T = smallestunsignedtype(m)
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
    VariedCostsMarket{T<:Unsigned}

Contains information about a college application market with varied application
costs.

# Constructors

```julia
VariedCostsMarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int}, H::Int)
```
Construct the `VariedCostsMarket` defined by admissions probabilities `f`,
utility values `t`, application costs `g`, and application limit `h`.

```jldoctest
julia> mkt = VariedCostsMarket([0.1, 0.5, 0.3, 0.1], [12, 3, 4, 13], [1, 2, 1, 1], 3)
VariedCostsMarket{UInt8}(0x04, [0.5, 0.3, 0.1, 0.1], [3, 4, 12, 13], [2, 1, 1, 1], 3, [1.5, 1.2, 1.2000000000000002, 1.3], [0.5, 0.7, 0.9, 0.9], UInt8[0x02, 0x03, 0x01, 0x04])
```

```julia
VariedCostsMarket(m)
```
Generate a random `VariedCostsMarket` with `m` schools.


# Internal API

This type contains the following fields:

 - `m`: Number of schools
 - `f`: Vector of admissions probabilities
 - `t`: Vector of utility values (must be sorted)
 - `g`: Vector of application costs
 - `H`: Budget to spend on applications
 - `ft = f .* t`
 - `omf = 1 .- f`
 - `perm`: How the input data were permuted to sort by `t`

`T` is the type of `h` and `m` and eltype of `perm`. Typically `T` is the
smallest `UInt` type for which `m < typemax(T)`. 

Internally, the schools must be indexed so that `t` is sorted ascending. `perm`
allows the original order to be recovered. See [`SameCostsMarket`](@ref).
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

    function VariedCostsMarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int}, H::Int)
        iscoherentmarket(f, t, g)
        T = smallestunsignedtype(length(f))
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
        T = smallestunsignedtype(m)
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

`SameCostsMarket` example (`h = 3` is irrelevant):

```jldoctest
julia> mkt = SameCostsMarket([0.1, 0.5, 0.3, 0.1], [12, 3, 4, 13], 3);

julia> round(valuation([1, 4], mkt), digits=2) # expected utility when applying to schools 1 and 4
2.38
```

`VariedCostsMarket` example (`g = [1, 2, 1, 1]` and `H = 4` are irrelevant):

```jldoctest
julia> mkt = VariedCostsMarket([0.1, 0.5, 0.3, 0.1], [12, 3, 4, 13], [1, 2, 1, 1], 4);

julia> round(valuation([1, 4], mkt), digits=2) # expected utility when applying to schools 1 and 4
2.38
```
"""
function valuation(
    X::AbstractVector{<:Integer},
    mkt::Union{SameCostsMarket{T},VariedCostsMarket{T}};
    invp::AbstractVector{T}=invperm(mkt.perm)
)::Float64 where {T<:Unsigned}
    isempty(X) && return 0.0

    X = T.(X)

    X[:] = sort(invp[X])

    h = T(length(X))

    if h > 1
        res = 0.0
        cp = reverse(cumprod(reverse(mkt.omf[X[2:end]])))

        for j in T(1):T(h - 1)
            res += mkt.ft[X[j]] * cp[j]
        end

        res += mkt.ft[X[end]]

        return res
    else
        return mkt.ft[X[1]]
    end
end
