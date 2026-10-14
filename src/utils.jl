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
julia> using SymmetricTightBinding: split_complex

julia> t = [im,0]
2-element Vector{Complex{Int64}}:
 0 + 1im
 0 + 0im

julia> T = split_complex(t)
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

julia> split_complex(t)
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
"""
    orbital_positions(br::BandRep{D})                            --> Vector{DirectPoint{D}}
    orbital_positions(cbr::CompositeBandRep{D})
    orbital_positions(atbm::AbstractTightBindingModel{D})
    orbital_positions(aptbm::AbstractParameterizedTightBindingModel{D})

Return a list of positions associated with the orbitals of a `BandRep` or a
`CompositeBandRep`, following the canonical orbital ordering of [`OrbitalOrdering`](@ref).
For a `CompositeBandRep`, the orbitals of each featured `BandRep` are concatenated, in the
ordering of their coefficients. For coefficients greater than 1, the positions are repeated
`cᵢ-1` times.
"""
function orbital_positions(br::BandRep{D}) where D
    # positions must be concrete, so we disallow free parameters (cf. `primitivized_orbit`)
    return [DirectPoint{D}(constant(o.wp)) for o in OrbitalOrdering(br; allow_free = false)]
end

function orbital_positions(cbr::CompositeBandRep{D}) where D
    N = occupation(cbr)
    positions = Vector{DirectPoint{D}}(undef, N)
    j = 0
    for (i, cᵢ) in enumerate(cbr.coefs)
        iszero(cᵢ) && continue
        positionsᵢ = orbital_positions(cbr.brs[i])
        Nᵢ = length(positionsᵢ)
        for _ in 1:Int(cᵢ) # add the positions once for each copy of the BR
            positions[j+1:j+Nᵢ] .= positionsᵢ
            j += Nᵢ
        end
    end
    j == N || error("inconsistent size calculation of `positions` vector")

    return positions
end

"""
    primitivized_orbit(br::BandRep{D}; allow_free::Bool = false)
                                                         --> Vector{WyckoffPosition{D}}

Return the orbit of the Wyckoff position associated with the band representation `br`, with
coordinates referred to the primitive basis. The order of the orbit is that of
`orbit(group(br))`.

The following checks are made, producing an error if violated, as the implementation
assumes:
1. There are no free parameters associated with the Wyckoff position; skipped if
   `allow_free = true` (free parameters are carried along symbolically by e.g., the
   enumeration of hopping orbits, but must be pinned to obtain concrete positions).
2. For every position, its coordinates, referred to the primitive basis, is in the range
   [0,1); i.e., every position lies in the parallelepiped primitive unit cell [0,1)ᴰ.
"""
function primitivized_orbit(br::BandRep{D}; allow_free::Bool = false) where D
    # we only want to include the Wyckoff positions in the primitive cell - but the default
    # listings from `spacegroup` include operations that are "centering translations";
    # fortunately, the orbit returned for a `BandRep` do not include these redundant
    # operations - but is still specified in a conventional basis. So, below, we change the
    # positions from a conventional to a primitive basis
    cntr = centering(num(br), D)
    wps = primitivize.(orbit(group(br)), cntr)
    for wp in wps
        if !allow_free && !iszero(free(wp))
            error(lazy"encountered Wyckoff position $wp with free parameters: not allowed")
        end
        if any(rᵢ -> rᵢ < 0 || rᵢ ≥ 1, constant(wp))
            error(lazy"encountered Wyckoff position $wp with primitive coordinates outside [0,1): this is inconsistent with implementation expectations, please file a bug report")
        end
    end
    return wps
end

"""
    primitivized_generators(br::BandRep{D}) --> Vector{<:AbstractOperation{D}}
    primitivized_generators(gens::AbstractVector{<:AbstractOperation{D}}, sgnum::Integer)

Return the generators of the space group of `br` in a *primitive* basis; or, alternatively,
the generators `gens`, given in the *conventional* basis of the space group `sgnum`.

Ē (the barred identity, i.e., a 2π rotation) is omitted from double group generators: it
imposes no constraint on a tight-binding model, since its representation is `-𝟙` (which
cancels in `ρ M ρ†`) and it leaves **k** invariant.

See also Crystalline.jl's `generators` function, which returns the generators in a
*conventional* basis.
"""
function primitivized_generators(br::BandRep{D}) where D
    sgnum = num(br)
    gens = generators(sgnum, isspinful(br) ? DSpaceGroup{D} : SpaceGroup{D})
    return primitivized_generators!(gens, sgnum)
end
function primitivized_generators(
    gens::AbstractVector{<:AbstractOperation{D}}, sgnum::Integer
) where D
    primitivized_generators!(copy(gens), sgnum)
end

"""
    primitivized_generators!(gens::AbstractVector{<:AbstractOperation{D}}, sgnum::Integer)

In-place version of [`primitivized_generators(gens, sgnum)`](@ref), mutating and returning
`gens`.
"""
function primitivized_generators!(
    gens::AbstractVector{O}, sgnum::Integer
) where {D, O<:AbstractOperation{D}}
    if isspinful(O)
        gens = filter!(g -> !(isbarred(g) && isone(SymOperation(g))), gens) # omit Ē
    end
    cntr = centering(sgnum, D)
    return cntr ∈ ('P', 'p') ? gens : map!(Base.Fix2(primitivize, cntr), gens)
end

# ---------------------------------------------------------------------------------------- #

"""
    pin_free!(
        brs::Collection{BandRep{D}},
        idx2αβγ::Pair{Int, <:AbstractVector{<:Real}}
    )

    pin_free!(
        brs::Collection{BandRep{D}},
        idx2αβγs::AbstractVector{<:Pair{Int, <:AbstractVector{<:Real}}}
    )

For `idx2αβγ = idx => αβγ`, update `brs[idx]` such that the free parameters of its
associated Wyckoff positions are pinned to `αβγ`.

A vector of pairs `idx2αβγs` can also be provided, to pin multiple distinct band
representations.

See also [`pin_free`](@ref) for non-mutated input.
"""
function pin_free!(
    brs::Collection{<:BandRep},
    idx2αβγs::AbstractVector{<:Pair{Int, <:AbstractVector{<:Real}}},
)
    foreach(Base.Fix1(pin_free!, brs), idx2αβγs)
    return brs
end

function pin_free!(
    brs::Collection{<:BandRep},
    idx2αβγ::Pair{Int, <:AbstractVector{<:Real}},
)
    idx, αβγ = idx2αβγ
    checkbounds(Bool, brs, idx) || error("index $idx out of bounds for `brs`")
    br = @inbounds brs[idx]
    @inbounds brs[idx] = pin_free(br, αβγ)
    return brs
end

"""
    pin_free(br::BandRep{D}, αβγ::AbstractVector{<:Real}) where D

Pin the free parameters of the Wyckoff position associated with the band representation `br`
to the values in `αβγ`.

Returns a new band representation with all other properties, apart from the Wyckoff
position, identical to (and sharing memory with) `br`.

Note that the associated orbit of the Wyckoff position will be automatically adjusted to
ensure that each position in the orbit lies within the primitive unit cell [0,1)ᴰ. That is,
if a choice of αβγ sends a position in the orbit outside the primitive unit cell, the
position will be adjusted by integer lattice translations to lie within.
"""
function pin_free(br::BandRep{D}, αβγ::AbstractVector{<:Real}) where D
    length(αβγ) == D || error(DimensionMismatch("length(αβγ) ≠ D"))
    if iszero(free(position(br)))
        error("attempting to pin a band representation without any free parameters")
    end

    wp = position(br)
    rv = parent(wp)
    rv_pin = RVec{D}(rv(αβγ))
    wp_pin = WyckoffPosition(wp.mult, wp.letter, rv_pin)

    siteir = br.siteir
    siteg = group(siteir)
    siteg_pin = typeof(siteg)(siteg.num, wp_pin, siteg.operations, siteg.cosets)

    # if the Wyckoff position that was picked is not in the primitive unit cell - or even if
    # a position in its orbit is not - we need to adjust the Wyckoff positions to lie inside
    # the primitive cell [0, 1)ᴰ: generally, this entails adjusting the choice of cosets
    # (which generate the orbit) and potentially also the site group operations for the
    # representative Wyckoff position
    orbit_pin = orbit(siteg_pin)
    cntr = centering(num(br), D)
    in_primitive_cell = all(orbit_pin) do rv_pin′
        rv_pin′_primitive = primitivize(rv_pin′, cntr)
        all(rᵢ -> 0 ≤ rᵢ < 1, constant(rv_pin′_primitive))
    end
    if !in_primitive_cell
        siteg_pin, _ = Crystalline.reduce_orbits_and_cosets(siteg_pin)
    end

    siteir_pin = typeof(siteir)(
        siteir.cdml,
        siteg_pin,
        siteir.matrices,
        siteir.reality,
        siteir.iscorep,
        siteir.pglabel,
    )

    return BandRep(siteir_pin, br.n, br.timereversal)
end

function reciprocal_translation_phase(
    positions::AbstractVector{DirectPoint{D}},
    k::ReciprocalPointLike{D},
) where D
    expsv = cispi.(dot.(Ref(-2 .* k), positions))
    return Diagonal(expsv)
end