function iscoherentmarket(f::Vector{Float64}, t::Vector{<:Real})
    length(f) == length(t) || throw(DimensionMismatch("`f` and `t` must have same length"))
    all(0 .< f .≤ 1) || throw(DomainError(f, "Must have `all(0 .< f .≤ 1)`"))
    all(0 .≤ t) || throw(DomainError(f, "Must have `all(0 .≤ t)`"))
end


function iscoherentmarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int})
    length(f) == length(t) || throw(DimensionMismatch("`f` and `t` must have same length"))
    length(f) == length(g) || throw(DimensionMismatch("`f` and `g` must have same length"))
    all(0 .< f .≤ 1) || throw(DomainError(f, "Must have `all(0 .< f .≤ 1)`"))
    all(0 .≤ t) || throw(DomainError(f, "Must have `all(0 .≤ t)`"))
    all(0 .≤ g) || throw(DomainError(f, "Must have `all(0 .≤ g)`"))
end


"""
Contains information about a college application market.
"""
abstract type Market end


"""
    SameCostsMarket{U<:Real}

Contains information about a college application market with identical
application costs.

# Constructors

```julia
SameCostsMarket(f::Vector{Float64}, t::Vector{U}, h::Integer) where {U<:Real}
```
Construct the `SameCostsMarket` defined by admissions probabilities `f`,
utility values `t`, and application limit `h`.

```julia-repl
julia> SameCostsMarket([0.3, 0.2, 0.05], [2, 3, 4], 2)
SameCostsMarket{Int64}(3, [0.3, 0.2, 0.05], [2, 3, 4], 2, [1, 2, 3])
```

```julia
SameCostsMarket(m)
```
Generate a random `SameCostsMarket{Int}` with `m` schools. `f` and `t` correlate negatively, mimicking
a realistic market.

# Internal API

This type contains the following fields:

 - `m`: Number of schools
 - `f`: Vector of admissions probabilities
 - `t`: Vector of utility values (must be sorted)
 - `h`: Number of schools student is allowed to apply to
 - `perm`: How the input data were permuted to sort by `t`

`U` is the eltype of `t`.

Internally, the schools must be indexed so that `t` is sorted ascending. `perm`
allows the original order to be recovered:

```julia-repl
julia> f = rand(3); t = rand(3); mkt = SameCostsMarket(f, t, 2);

julia> t[mkt.perm] == mkt.t
true

julia> t == mkt.t[invperm(mkt.perm)]
true
```
"""
mutable struct SameCostsMarket{U<:Real} <: Market
    m::Int
    f::Vector{Float64}
    t::Vector{U}
    h::Int
    perm::Vector{Int}

    function SameCostsMarket(f::Vector{Float64}, t::Vector{U}, h::Integer) where {U<:Real}
        iscoherentmarket(f, t)
        perm = sortperm(t)

        return new{U}(length(f), f[perm], t[perm], min(h, length(f)), perm)
    end

    function SameCostsMarket(m::Integer)
        t = ceil.(Int, -10 * log.(rand(m))) |> sort
        f = inv.(t .+ 10 * rand(m))
        perm = 1:m
        return new{Int}(m, f[perm], t[perm], m ÷ 2, perm)
    end
end


"""
    VariedCostsMarket

Contains information about a college application market with varied application
costs.

# Constructors

```julia
VariedCostsMarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int}, H::Int)
```
Construct the `VariedCostsMarket` defined by admissions probabilities `f`,
utility values `t`, application costs `g`, and application limit `h`.

```julia-repl
julia> mkt = VariedCostsMarket([0.1, 0.5, 0.3, 0.1], [12, 3, 4, 13], [1, 2, 1, 1], 3)
VariedCostsMarket(4, [0.5, 0.3, 0.1, 0.1], [3, 4, 12, 13], [2, 1, 1, 1], 3, [2, 3, 1, 4])
```

```julia
VariedCostsMarket(m)
```
Generate a random `SameCostsMarket{Int}` with `m` schools. `f` and `t` correlate negatively, mimicking
a realistic market. Entries of `g` are random integers in `5:10`.

# Internal API

This type contains the following fields:

 - `m`: Number of schools
 - `f`: Vector of admissions probabilities
 - `t`: Vector of utility values (must be sorted)
 - `g`: Vector of application costs
 - `H`: Budget to spend on applications
 - `perm`: How the input data were permuted to sort by `t`

Internally, the schools must be indexed so that `t` is sorted ascending. `perm`
allows the original order to be recovered:

```julia-repl
julia> f = rand(3); t = rand(3); g = rand(5:10, 3); mkt = VariedCostsMarket(f, t, g, 10);

julia> t[mkt.perm] == mkt.t
true

julia> t == mkt.t[invperm(mkt.perm)]
true
```
"""
mutable struct VariedCostsMarket <: Market
    m::Int
    f::Vector{Float64}
    t::Vector{Int}
    g::Vector{Int}
    H::Int
    perm::Vector{Int}

    function VariedCostsMarket(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int}, H::Int)
        iscoherentmarket(f, t, g)
        perm = sortperm(t)

        return new(length(f), f[perm], t[perm], g[perm], min(H, sum(g)), perm)
    end

    function VariedCostsMarket(m::Integer)
        t = ceil.(Int, -10 * log.(rand(m))) |> sort
        f = inv.(t .+ 10 * rand(m))
        g = rand(5:10, m)
        H = sum(g) ÷ 2
        perm = 1:m

        return new(m, f[perm], t[perm], g[perm], H, perm)
    end
end


"""
    Market(f, t, h)

Equivalent to [`SameCostsMarket(f, t, h)`](@ref).
"""
function Market(f::Vector{Float64}, t::Vector{<:Real}, h::Int)
    return SameCostsMarket(f, t, h)
end


"""
    Market(f, t, g, H)

Equivalent to [`VariedCostsMarket(f, t, g, H)`](@ref).
"""
function Market(f::Vector{Float64}, t::Vector{Int}, g::Vector{Int}, H::Int)
    return VariedCostsMarket(f, t, g, H)
end


function valuation_nopermute_sorted(X::AbstractVector{<:Integer}, mkt::Market;)::Float64
    isempty(X) && return 0.0
    # @assert issorted(X)
    res = 0.0

    for j in X
        res = (1 - mkt.f[j]) * res + mkt.f[j] * mkt.t[j]
    end

    return res
end


"""
    valuation(X, mkt)

Return the valuation of the portfolio `X` for the market `mkt`, which may be either 
a `SameCostsMarket` or a `VariedCostsMarket`.

`SameCostsMarket` example (`h = 3` is irrelevant):

```julia-repl
julia> mkt = SameCostsMarket([0.1, 0.5, 0.3, 0.1], [12, 3, 4, 13], 3);

julia> round(valuation([1, 4], mkt), digits=2) # expected utility when applying to schools 1 and 4
2.38
```

`VariedCostsMarket` example (`g = [1, 2, 1, 1]` and `H = 4` are irrelevant):

```julia-repl
julia> mkt = VariedCostsMarket([0.1, 0.5, 0.3, 0.1], [12, 3, 4, 13], [1, 2, 1, 1], 4);

julia> round(valuation([1, 4], mkt), digits=2) # expected utility when applying to schools 1 and 4
2.38
```
"""
function valuation(
    X::AbstractVector{<:Integer},
    mkt::Market;
    invp::AbstractVector{<:Integer} = invperm(mkt.perm),
)::Float64
    return valuation_nopermute_sorted(sort(invp[X]), mkt)
end
