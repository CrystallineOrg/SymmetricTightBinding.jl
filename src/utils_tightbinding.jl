"""
    split_complex(t::Vector{<:Number}) -> Matrix{Real}

Consider `αt` where `α ∈ ℂ` and `t ∈ ℂⁿ` and build from `t` a matrix representation
`T` that allows access to the real and imaginary parts of the product `αt` without using
complex numbers by splitting α into a real 2-vector of its real and imaginary parts.

In particular, let ``α = αᴿ + iαᴵ`` and ``t = tᴿ + itᴵ`` with `αᴿ, αᴵ ∈ ℝ` and
``tᴿ, tᴵ ∈ ℝⁿ``, then ``αt`` can be rewritten as

```math
αt = (αᴿ + iαᴵ)(tᴿ + itᴵ)
   = (αᴿtᴿ - αᴵtᴵ) + i(αᴿtᴵ + αᴵtᴿ)
   = [tᴿ, tᴵ]ᵀ [αᴿ, αᴵ] + i [tᴵ, tᴿ]ᵀ [αᴿ, αᴵ]
```

Then, defining `T = [tᴿ -tᴵ; tᴵ tᴿ]`, the above product can then be reexpressed as:
``Re(αt) = αᴿtᴿ - αᴵtᴵ =`` `(T * [αᴿ; αᴵ])[1:n]` and ``Im(αt) = αᴿtᴵ + αᴵtᴿ =``
`(T * [αᴿ; αᴵ])[n+1:2n]`.
I.e., the "upper half" of the product `T * [real(α), imag(α)]` is `real(α * t)` and the 
"lower half" is `imag(αt)`.

This functionality is used to avoid complex numbers in amplitude basis coefficients, which
simplifies the application of time-reversal symmetry and hermiticity.

## Examples

```julia
julia> using TETB: split_complex

julia> t = [im,0]
2-element Vector{Complex{Int64}}:
 0 + 1im
 0 + 0im

julia> T = TETB.split_complex(t)
4×2 Matrix{Int64}:
 0  -1
 0   0
 1   0
 0   0

julia> α = 0.5+0.2im; αv = [real(α), imag(α)];

julia> (T * αv)[1:2] == real(α*t) && (T * αv)[3:4] == imag(α*t)
```

```julia
julia> t = [1,im]
2-element Vector{Complex{Int64}}:
 1 + 0im
 0 + 1im

julia> TETB.split_complex(t)
4×2 Matrix{Int64}:
 1   0
 0  -1
 0   1
 1   0
```
"""
function split_complex(t::AbstractVector{<:Number})
    re_t, im_t = reim(t)
    return [re_t -im_t; im_t re_t] # == [real(t) real(im * t); imag(t) imag(im * t)]
end

"""
    inversion(::Val{D}) --> SymOperation{D}

Return the inversion operation in dimension `D`.
"""
inversion(::Val{3}) = S"-x,-y,-z"
inversion(::Val{2}) = S"-x,-y"
inversion(::Val{1}) = S"-x"
inversion(::Val) = error("unsupported dimension")

## --------------------------------------------------------------------------------------- #
# Extracting a list of positions associated with our convention for orbital ordering of a
# `NewBandRep` or a `CompositeBandRep`. For a `NewBandRep`, the orbitals are arranged such
# that the first `irdim(br.siteir)` orbitals associate to the first element of the orbit
# of its Wyckoff positions; the next `irdim(br.siteir)` orbitals associate to the second
# element of the orbit, and so on. For a `CompositeBandRep`, the orbitals of each
# `NewBandRep` are concatenated, in the order of their coefficients. For coefficients
# greater than 1, the positions are repeated `cᵢ-1` times.

function orbital_positions(br::NewBandRep{D}) where D
    dim = irdim(br.siteir)
    wps = primitivized_orbit(br)
    positions = Vector{DirectPoint{D}}(undef, length(wps) * dim)
    for (m, wp) in enumerate(wps)
        for j in ((m-1)*dim+1):(m*dim)
            positions[j] = wp
        end
    end
    return positions
end

function orbital_positions(cbr::CompositeBandRep{D}) where D
    N = occupation(cbr)
    positions = Vector{DirectPoint{D}}(undef, N)
    j = 0
    for (i, cᵢ) in enumerate(cbr.coefs)
        iszero(cᵢ) && continue
        br = cbr.brs[i]

        dim = irdim(br.siteir)
        wps = primitivized_orbit(br)
        Nᵢ = length(wps) * dim
        for (m, wp) in enumerate(wps)
            for j′ in ((m-1)*dim+1):(m*dim)
                positions[j+j′] = wp
            end
        end
        j += Nᵢ

        # if cᵢ > 1, we add the just-added positions `cᵢ-1` times more
        for _ in 1:(Int(cᵢ)-1)
            @views positions[j+1:j+Nᵢ] .= positions[j-Nᵢ+1:j]
            j += Nᵢ
        end
    end
    j == N || error("inconsistent size calculation of `positions` vector")

    return positions
end

"""
    primitivized_orbit(br::NewBandRep{D}) where D

Return the orbit of the Wyckoff position associated with the band representation `br`.
The coordinates of positions in the orbit are given relative to the primitive unit cell.

Positions are returned as a `Vector{DirectPoint{D}}`.

The following checks are made, producing an error if violated:
1. There are no free parameters associated with the Wyckoff position.
2. For every position, its coordinates, referred to the primitive basis, is in the range
   [0,1); i.e., every position lies in the paralleliped primitive unit cell [0,1)ᴰ.
"""
function primitivized_orbit(br::NewBandRep{D}) where D
    wps = orbit(group(br))
    cntr = centering(num(br), D)
    wps′_pts = Vector{DirectPoint{D}}(undef, length(wps))
    for (m, wp) in enumerate(wps)
        wp′ = primitivize(wp, cntr)
        if !iszero(free(wp′))
            error(lazy"encountered Wyckoff position $wp with free parameters: not allowed")
        end
        wp′_cnst = constant(wp′)
        if any(rᵢ -> rᵢ < 0 || rᵢ ≥ 1, wp′_cnst)
            error(
                lazy"encountered Wyckoff position $wp (conventional coordinates) with primitive coordinates $wp′_cnst outside [0,1): this inconsistent with implementation expectations, please file a bug report",
            )
        end
        wps′_pts[m] = DirectPoint{D}(wp′_cnst)
    end
    return wps′_pts
end
