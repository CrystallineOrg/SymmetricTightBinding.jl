# Composite tight-binding models for generally non-Hermitian Hamiltonians, built from a
# Hermitian and an anti-Hermitian `TightBindingModel`. Split out from `types.jl` since the
# content is self-contained; the associated `show` methods live in `show.jl`.

# ---------------------------------------------------------------------------------------- #

const HorAH_TBT{D} = Union{TightBindingTerm{D, HERMITIAN}, TightBindingTerm{D, ANTIHERMITIAN}}

"""
    CompositeTightBindingModel{D} <: AbstractTightBindingModel{HorAH_TBT{D}}

Composite tight-binding model for generally non-Hermitian Hamiltonians, featuring both
Hermitian and anti-Hermitian terms, stored as separate `TightBindingModel`s.

A `CompositeTightBindingModel` is created from a Hermitian and an anti-Hermitian
`TightBindingModel`, either via `CompositeTightBindingModel(tbm_h, tbm_a)` or via the
shorthand `tbm_h + tbm_a` (in both cases, the argument order is immaterial). Both models
must be founded on the same composite band representation, e.g.:

```jl
tbm_h = tb_hamiltonian(cbr, Rs, Val(HERMITIAN))
tbm_a = tb_hamiltonian(cbr, Rs, Val(ANTIHERMITIAN))
ctbm = tbm_h + tbm_a
```

Iterating a `CompositeTightBindingModel` traverses the Hermitian terms first, then the
anti-Hermitian terms; the element type is correspondingly the 2-union of
`TightBindingTerm{D, S}` spanning `S = HERMITIAN` and `ANTIHERMITIAN`.

To associate a set of coefficients to each term, see
[`ParameterizedCompositeTightBindingModel`](@ref).

## Fields
- `tbm_h :: TightBindingModel{D, HERMITIAN}`: the Hermitian part of the model
- `tbm_a :: TightBindingModel{D, ANTIHERMITIAN}`: the anti-Hermitian part of the model
"""
struct CompositeTightBindingModel{D} <: AbstractTightBindingModel{HorAH_TBT{D}}
    tbm_h::TightBindingModel{D, HERMITIAN}
    tbm_a::TightBindingModel{D, ANTIHERMITIAN}
    function CompositeTightBindingModel{D}(
        tbm_h::TightBindingModel{D, HERMITIAN},
        tbm_a::TightBindingModel{D, ANTIHERMITIAN},
    ) where D
        # validate inputs
        tbm_h.N == tbm_a.N || error("input models must have the same number of orbitals")
        if CompositeBandRep(tbm_h) ≠ CompositeBandRep(tbm_a)
            error("input models must be based on the same composite band representation")
        end
        if orbital_positions(tbm_h) ≠ orbital_positions(tbm_a)
            error("input models must have the same orbital positions")
        end
        return new{D}(tbm_h, tbm_a)
    end
end
function CompositeTightBindingModel(
    tbm_h::TightBindingModel{D, HERMITIAN},
    tbm_a::TightBindingModel{D, ANTIHERMITIAN},
) where D
    return CompositeTightBindingModel{D}(tbm_h, tbm_a)
end
function Base.:+(
    tbm_h::TightBindingModel{D, HERMITIAN},
    tbm_a::TightBindingModel{D, ANTIHERMITIAN}
) where D
    return CompositeTightBindingModel{D}(tbm_h, tbm_a)
end

# inverted argument order variants
# NB: the `{D}` entry below must be spliced as an _expression_ (`:(…{D})`), not as a
#     `Symbol("…{D}")`: the latter defines a stray function literally named `var"…{D}"`
#     rather than a parametric constructor method
for f in (:(CompositeTightBindingModel{D}), :CompositeTightBindingModel, :(Base.:+))
    @eval function $f(
        tbm_a::TightBindingModel{D, ANTIHERMITIAN},
        tbm_h::TightBindingModel{D, HERMITIAN},
    ) where D
        return $f(tbm_h, tbm_a)
    end
end

# AbstractTightBindingModel interface
hermiticity(::CompositeTightBindingModel) = NONHERMITIAN
orbital_positions(ctbm::CompositeTightBindingModel) = orbital_positions(ctbm.tbm_h)
Crystalline.CompositeBandRep(ctbm::CompositeTightBindingModel) = CompositeBandRep(ctbm.tbm_h)
orbital_count(ctbm::CompositeTightBindingModel) = orbital_count(ctbm.tbm_h)
Crystalline.dim(::CompositeTightBindingModel{D}) where D = D

# AbstractVector interface
# NB: `eltype` follows from the `AbstractVector{HorAH_TBT{D}}` supertype; no method needed.
# NB: `similar` is intentionally _not_ extended: it is handed only the requested `dims`,
#     never the indices, so it cannot know how a request splits across the Hermitian and
#     anti-Hermitian parts. Index-preserving slicing is instead provided directly by the
#     `getindex(::CompositeTightBindingModel, ::AbstractVector)` methods below. Operations
#     that are genuinely not closed over the type — notably `vcat`, whose result interleaves
#     the two parts — consequently fall back to a plain `Vector{HorAH_TBT{D}}`.
Base.size(ctbm::CompositeTightBindingModel) = (length(ctbm.tbm_h)+length(ctbm.tbm_a),)
function Base.getindex(ctbm::CompositeTightBindingModel, i::Int)
    if 1 ≤ i ≤ length(ctbm.tbm_h)
        return @inbounds ctbm.tbm_h[i]
    elseif length(ctbm.tbm_h) < i ≤ length(ctbm)
        return @inbounds ctbm.tbm_a[i - length(ctbm.tbm_h)]
    else
        throw(BoundsError(ctbm, i))
    end
end
function Base.setindex!(ctbm::CompositeTightBindingModel, v, i::Int)
    if 1 ≤ i ≤ length(ctbm.tbm_h)
        return @inbounds setindex!(ctbm.tbm_h, v, i)
    elseif length(ctbm.tbm_h) < i ≤ length(ctbm)
        return @inbounds setindex!(ctbm.tbm_a, v, i - length(ctbm.tbm_h))
    else
        throw(BoundsError(ctbm, i))
    end
end
Base.IndexStyle(::Type{<:CompositeTightBindingModel}) = IndexLinear()

# Indexing by a vector of indices returns a `CompositeTightBindingModel` again, which is
# convenient for dropping terms from a model (e.g., `ctbm[[1,2,5]]` or `ctbm[setdiff(
# 1:length(ctbm), 3)]`). This is possible because the terms of `ctbm` are ordered with all
# Hermitian terms first, then all anti-Hermitian terms: any index set therefore partitions
# cleanly across the two parts.
# NB: `idxs` is required to be strictly increasing. A composite model can only ever store
#     its terms in the above canonical order, so a permuted `idxs` (e.g. `ctbm[[9,1]]`)
#     could not be honoured and would silently return a differently-ordered model than
#     requested; we error instead. Use `collect(ctbm)[idxs]` to obtain an arbitrarily
#     ordered `Vector` of terms.
function Base.getindex(
    ctbm::CompositeTightBindingModel{D},
    idxs::AbstractVector{<:Integer}
) where D
    for i in idxs
        checkbounds(ctbm, i)
    end
    if !issorted(idxs; lt = ≤) # i.e., not strictly increasing (also rejects duplicates)
        error("indexing a `CompositeTightBindingModel` requires strictly increasing \
               indices, since its terms are necessarily stored with all Hermitian terms \
               before all anti-Hermitian terms; use `collect(ctbm)[idxs]` to obtain an \
               arbitrarily ordered `Vector` of terms instead")
    end
    Nʰ = length(ctbm.tbm_h)
    idxsʰ = [i for i in idxs if i ≤ Nʰ]
    idxsᵃ = [i - Nʰ for i in idxs if i > Nʰ]
    return CompositeTightBindingModel{D}(ctbm.tbm_h[idxsʰ], ctbm.tbm_a[idxsᵃ])
end
# logical indexing (`ctbm[mask]`): defined separately since `Bool <: Integer`, which would
# otherwise make the method above misinterpret a mask as a vector of indices
function Base.getindex(ctbm::CompositeTightBindingModel, mask::AbstractVector{Bool})
    length(mask) == length(ctbm) || throw(BoundsError(ctbm, mask))
    return ctbm[findall(mask)]
end

function (ctbm::CompositeTightBindingModel{D})(cs::AbstractVector{<:Real}) where {D}
    return ParameterizedCompositeTightBindingModel{D}(ctbm, cs)
end

function (ctbm::CompositeTightBindingModel{D})(
    cs_h::AbstractVector{<:Real},
    cs_a::AbstractVector{<:Real}
) where {D}
    length(cs_h) == length(ctbm.tbm_h) || error("mismatched number of coefficients and terms for Hermitian part")
    length(cs_a) == length(ctbm.tbm_a) || error("mismatched number of coefficients and terms for anti-Hermitian part")
    return ParameterizedCompositeTightBindingModel{D}(ctbm, vcat(cs_h, cs_a))
end

## --------------------------------------------------------------------------------------- #

"""
    ParameterizedCompositeTightBindingModel{D} <: AbstractParameterizedTightBindingModel{D}

A coefficient-parameterized [`CompositeTightBindingModel`](@ref), that can be used as a
functor for evaluation at input momenta `k`.

A `ParameterizedCompositeTightBindingModel` is usually obtained by calling a
`CompositeTightBindingModel` `ctbm` as a functor, either with a single vector of
coefficients ordered as the terms of `ctbm` itself, i.e., the Hermitian coefficients
followed by the anti-Hermitian ones (`ctbm(cs)`), or with the two sets of coefficients
provided separately (`ctbm(cs_h, cs_a)`).

## Fields
- `tbm :: CompositeTightBindingModel{D}`: the underlying composite tight-binding model
- `cs :: Vector{Float64}`: coefficients of each term of `tbm`, in the same order; i.e.,
  `cs[1:length(tbm.tbm_h)]` parameterize the Hermitian terms and `cs[length(tbm.tbm_h)+1:
  end]` the anti-Hermitian terms
- `scratch :: Matrix{ComplexF64}`: scratch space for evaluation

## Evaluation
A `ParameterizedCompositeTightBindingModel`, `pctbm`, can be evaluated at any
`D`-dimensional momentum `k` by using `pctbm` as a functor, i.e., `pctbm(k)` returns the 
numerical (generally non-Hermitian) Hamiltonian matrix at momentum `k`.

!!! warning
    The returned matrix aliases the internal `scratch` buffer of `pctbm` and is overwritten
    by subsequent evaluations: `copy` it if it must outlive the next call.
"""
struct ParameterizedCompositeTightBindingModel{D} <: AbstractParameterizedTightBindingModel{D}
    tbm :: CompositeTightBindingModel{D}
    cs :: Vector{Float64} # coefficients of the tight-binding model
    scratch :: Matrix{ComplexF64} # scratch space for evaluation
    function ParameterizedCompositeTightBindingModel{D}(
        tbm :: CompositeTightBindingModel{D},
        cs :: AbstractVector{<:Real},
        scratch :: Matrix{ComplexF64} = Matrix{ComplexF64}(
                                            undef, orbital_count(tbm), orbital_count(tbm)),
    ) where D
        length(tbm) ≠ length(cs) && _throw_term_coef_length_mismatch(tbm, cs)
        N = orbital_count(tbm)
        size(scratch) ≠ (N, N) && _throw_scratch_size_mismatch(scratch, N)
        return new{D}(tbm, convert(Vector{Float64}, cs), scratch)
    end
end

hermiticity(::ParameterizedCompositeTightBindingModel) = NONHERMITIAN
function orbital_positions(pctbm::ParameterizedCompositeTightBindingModel)
    return orbital_positions(pctbm.tbm)
end
function Crystalline.CompositeBandRep(pctbm::ParameterizedCompositeTightBindingModel)
    return CompositeBandRep(pctbm.tbm)
end

function (pctbm::ParameterizedCompositeTightBindingModel{D})(
    k::ReciprocalPointLike{D},
    scratch::Matrix{ComplexF64} = pctbm.scratch,
) where {D}
    if length(k) ≠ D
        error("momentum `k` must be a $D-dimensional vector to match the model dimension")
    end
    tbm_h = pctbm.tbm.tbm_h
    tbm_a = pctbm.tbm.tbm_a
    N = orbital_count(pctbm)
    size(scratch) ≠ (N, N) && _throw_scratch_size_mismatch(scratch, N)

    H = scratch # grab & reset scratch space for evaluating Hamiltonian matrix
    fill!(H, 0.0)

    # evaluate each block of the Hamiltonian terms, multiply by coefficients, & store in `H`
    # NB: `cs[1:length(tbm_h)]` are the coefficients for the Hermitian terms (`tbm_h`), and
    #     `cs[length(tbm_h)+1:end]` are for the anti-Hermitian terms (`tbm_a`); below, we
    #     split the loop explicitly to restore type-stability
    for (tbt, c) in zip(tbm_h.terms, @view pctbm.cs[1:length(tbm_h)])
        evaluate_tight_binding_term!(tbt::TightBindingTerm{D, HERMITIAN}, k, c, H) # modifies `H` in-place
    end
    cs_offset = length(tbm_h)
    for (tbt, c) in zip(tbm_a.terms, @view pctbm.cs[cs_offset+1:end])
        evaluate_tight_binding_term!(tbt::TightBindingTerm{D, ANTIHERMITIAN}, k, c, H) # modifies `H` in-place
    end
    
    return H
end
