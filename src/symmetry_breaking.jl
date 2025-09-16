"""
    subduced_complement(tbm::TightBindingModel{D}, sgnumᴴ::Int; timereversal)
                                                        --> TightBindingModel{D}

Given a model `tbm` associated with a space group ``G``, determine the new, independent
tight-binding terms (i.e., the the orthogonal complement of terms) that become 
symmetry-allowed when the model's space group is reduced to a subgroup ``H ≤ G`` with space
group number `sgnumᴴ` and time-reversal symmetry `timereversal`.

Practically, the function answers the question: which new tight-binding terms become allowed
if the symmetry of the model is reduced from space group ``G`` to subgroup ``H``?

## Implementation

The function computes a basis of allowed tight-binding terms in the subgroup setting ``H``
by simply restricting the constraints in ``G`` to generators in ``H``. This gives a basis
for the tight-binding terms in the subduced ``G ↓ H`` setting. 
The space spanned by this basis is compared to the space spanned in the original model; in
particular new terms are identified as the orthogonal complement of the spaces associated
with ``G ↓ H`` relative to ``G``.

## Keywords
- `timereversal::Bool`: Specifies whether time-reversal symmetry is present in the
  subgroup ``H``. By default, the presence or absence is inherited from the original model
  `tbm`. Note that `timereversal` must be "lower or equal to" the time-reversal of the
  original model.

## Example

It is well-known that the Dirac point of graphene is gapped under mirror and time-reversal
symmetry breaking. We can see this by constructing a tight-binding model first for a model
of graphene (plane group ⋕17) and then subducing it to a setting without mirror symmetry
(plane group ⋕16) and without time-reversal symmetry (`timereversal = false`). First, we
construct the tight-binding model for graphene (via the (2a|A₁) band representation):
```julia-repl
julia> using SymmetricTightBinding, Crystalline

julia> brs = calc_bandreps(17, Val(2); timereversal = true);

julia> cbr = @composite brs[5]

julia> tbm = tb_hamiltonian(cbr, [[0,0], [1,0]])
```
Each of the 4 terms in this model is proportional to an identity matrix at K = (1/3, 1/3).
Using `subduced_complement`, we can find the new terms that appear if we imagine lowering
the symmetry from plane group ⋕17 to ⋕16 (which has no mirror symmetry) while also removing
time-reversal symmetry.
```julia-repl
julia> Δtbm = subduced_complement(tbm, 16; timereversal = false)
2-term 2×2 TightBindingModel{2} over (2b|A₁):
┌─
1. ⎡ i𝕖(δ₁)+i𝕖(δ₂)+i𝕖(δ₃)-i𝕖(δ₄)-i𝕖(δ₅)-i𝕖(δ₆)  0                                          ⎤
│  ⎣ 0                                          -i𝕖(δ₁)-i𝕖(δ₂)-i𝕖(δ₃)+i𝕖(δ₄)+i𝕖(δ₅)+i𝕖(δ₆) ⎦
└─ (2b|A₁) self-term:  δ₁=[1,0], δ₂=[0,1], δ₃=[-1,-1], δ₄=-δ₁, δ₅=-δ₂, δ₆=-δ₃
┌─
2. ⎡ 0                                       𝕖(δ₁)+𝕖(δ₂)+𝕖(δ₃)-𝕖(δ₇)-𝕖(δ₈)-𝕖(δ₉) ⎤
│  ⎣ 𝕖(δ₄)+𝕖(δ₅)+𝕖(δ₆)-𝕖(δ₁₀)-𝕖(δ₁₁)-𝕖(δ₁₂)  0                                   ⎦
└─ (2b|A₁) self-term:  δ₁=[4/3,-1/3], δ₂=[1/3,5/3], δ₃=[-5/3,-4/3], δ₄=-δ₁, δ₅=-δ₂, δ₆=-δ₃, δ₇=[1/3,-4/3], δ₈=[-5/3,-1/3], δ₉=[4/3,5/3], δ₁₀=-δ₇, δ₁₁=-δ₈, δ₁₂=-δ₉
```
The first of the of these terms is not diagonal at K and so opens a gap at the Dirac point:
```julia-repl
julia> Δtbm[1](ReciprocalPoint(1/3, 1/3))
2×2 Matrix{ComplexF64}:
 5.19615+1.73195e-14im       0.0+0.0im
     0.0+0.0im          -5.19615+1.43774e-14im
```

### Adding symmetry-breaking terms to the original model
To build a "complete" model, with both the original and symmetry-breaking terms, use `vcat`:
```julia-repl
julia> tbm′ = vcat(tbm, Δtbm); length(tbm′) == length(tbm) + length(Δtbm)
true
```

## Limitations
The subgroup ``H`` must be a volume-preserving subgroup of the original group ``G``. I.e.
``H`` must be a translationen-gleiche subgroup of ``G`` (or ``G`` itself), and there must
exist a transformation from ``G`` to ``H`` that preserves volume (i.e., has
`det(t.P) == 1` for `t` denoting an element returned by Crystalline.jl's
`conjugacy_relations`).
"""
function subduced_complement(tbm::TightBindingModel{D}, sgnumᴴ::Int; kws...) where D
    sgnumᴳ = num(tbm.cbr)
    gr = maximal_subgroups(sgnumᴳ, SpaceGroup{D})
    ts = conjugacy_relations(gr, sgnumᴳ, sgnumᴴ)
    # note: it doesn't matter which of the conjugacy transforms we pick - we just need to
    # be able to transform the generators of H to the setting of G - in the end, we want
    # to present the results in our original setting (G), so it doesn't matter _which_ H
    # setting we imagine starting from (so long that it preserves volume)
    i = findfirst(ts) do t′
        det(t′.P) == 1
    end
    if isnothing(i)
        error("could not find a volume-preserving transformation to subgroup: ensure \
               that the groups have identical centerings")
    end
    t = ts[something(i)]  # note: `t` = (P|p) will take G to H - we want the opposite
    Pᴴ²ᴳ = inv(t.P)       # opposite transform: (Pᴴ²ᴳ|pᴴ²ᴳ) = (P|p)⁻¹
    pᴴ²ᴳ = -Pᴴ²ᴳ*t.p

    _gensᴴ = generators(sgnumᴴ, SpaceGroup{D}) # in H setting
    gensᴴ = transform.(_gensᴴ, Ref(Pᴴ²ᴳ), Ref(pᴴ²ᴳ))

    return subduced_complement(tbm, gensᴴ; kws...)
end

function subduced_complement(
    tbm::TightBindingModel{D},
    gensᴴ::AbstractVector{SymOperation{D}};
    timereversal::Bool = first(tbm.cbr.brs).timereversal, # ← whether H has time-reversal
) where D
    timereversalᴳ = first(tbm.cbr.brs).timereversal
    if timereversalᴳ == false && timereversal == true
        error("requested subgroup `timereversal = false`, but original model was built with time-reversal present; input for H must maintain or reduce symmetry")
    end

    # we need to go through the terms of `tbm` in "groups of the same orbit" - each orbit
    # will have some coefficient basis, and it is this basis we need to compare. So first,
    # we figure out the groupings into these orbits by just looking at `tbt.block.h_orbit`
    grouped_orbits_idxs = UnitRange{Int}[]
    current_h_orbit = first(tbm.terms).block.h_orbit
    i₁ = 1; i₂ = 0
    for tbt in tbm.terms
        if current_h_orbit == tbt.block.h_orbit
            i₂ += 1
        else
            push!(grouped_orbits_idxs, i₁:i₂)
            current_h_orbit = tbt.block.h_orbit
            i₁ = i₂ = i₂+1
        end
    end
    i₂ == length(tbm.terms) && push!(grouped_orbits_idxs, i₁:i₂)

    # now we can compute a new coefficient basis in H and compare with our original basis,
    # progressing "group by group"
    complement_tbs = TightBindingTerm{D}[]
    for idxs in grouped_orbits_idxs
        tbt = tbm.terms[first(idxs)]
        tbb = tbt.block
        # first, compute basis of coefficients for new subset of generators (`gensᴴ`)
        tₐᵦ_basis_reimᴴ_vs = _obtain_basis_free_parameters(
            tbb.h_orbit,
            tbb.br1,
            tbb.br2,
            tbb.ordering1, 
            tbb.ordering2, 
            tbb.Mm,
            gensᴴ,
            timereversal,
            tbt.block_ij[1] == tbt.block_ij[2], #= .diagonal_block =# 
            tbt.hermiticity == ANTIHERMITIAN,   #= .antihermitian =# 
        )
        # check output dimensions
        if length(tₐᵦ_basis_reimᴴ_vs) < length(idxs)
            error("unexpectedly found lower-dimensional basis space for model in subduced \
                   group; unexpected and unhandled - make sure the generators are a subset \
                   of the original generators (i.e., that fewer constraints apply than \
                   originally)")
        elseif length(tₐᵦ_basis_reimᴴ_vs) == length(idxs)
            continue # basis must then be unchanged; nothing to add for this index group
        end

        # get "original" coefficient basis in G from `tbm[idxs]
        tₐᵦ_basis_reimᴴ = stack(tₐᵦ_basis_reimᴴ_vs)
        tₐᵦ_basis_reimᴳ = Matrix{Float64}(undef, length(tbb.t), length(idxs))
        for (n, i) in enumerate(idxs)
            tbbᵢ = tbm.terms[i].block
            tₐᵦ_basis_reimᴳ[:,n] .= tbbᵢ.t
        end

        # find the orthogonal complement of the span of `tₐᵦ_basis_reimᴴ` relative to
        # `tₐᵦ_basis_reimᴳ` using QR factorization (i.e., find a basis for the space that is
        # in H but not in G)
        Qᴳ = Matrix(qr(tₐᵦ_basis_reimᴳ).Q) # columns of Qᴳ form orthonormal basis for G's coefs
        Pᴳ = Qᴳ * transpose(Qᴳ) # projection onto G space
        Pᵪᴳ = I - Pᴳ            # projection onto orthogonal complement of G space
        tₐᵦ_basis_reim_ᴴᵪᴳ = Pᵪᴳ * tₐᵦ_basis_reimᴴ # ortho. complement of H coefs rel. to G

        # now, extract a basis for the span of `tₐᵦ_basis_reim_ᴴᵪᴳ` (vectors could be near
        # zero or simply linearly dependent): rather than using QR, we use the SVD, so we
        # can avoid picking up basis elements that are just due to accumulated (floating
        # point, e.g.) errors
        Uᴴᵪᴳ, σs, _ = svd(tₐᵦ_basis_reim_ᴴᵪᴳ) # = U*Σ*Vᵀ w/ Σ = Diagonal(σs)
        Nᴴ = size(tₐᵦ_basis_reimᴴ, 2) - size(tₐᵦ_basis_reimᴳ, 2)
        tₐᵦ_basis_reim_ᴴᵪᴳ′ = Matrix{Float64}(undef, length(tbb.t), Nᴴ)
        for (n, (u, σ)) in enumerate(zip(eachcol(Uᴴᵪᴳ), σs))
            n > Nᴴ && continue # there should be exactly Nᴴ non-zero singular values
            if σ < NULLSPACE_ATOL_DEFAULT
                # make sure we don't have any surprises & verify our assumption above
                error(LazyString("unexpectedly found near-zero singular value for a SVD \
                      column vector that ought to have been a proper basis vector (Nᴴ = ",
                      Nᴴ, " σs = ", σs, ")"))
            end
            # keep this vector
            tₐᵦ_basis_reim_ᴴᵪᴳ′[:,n] .= u
        end

        # finally, make the basis vectors we have now "pretty"
        tₐᵦ_basis_reim_ᴴᵪᴳ′_sparsified = _poormans_sparsification(tₐᵦ_basis_reim_ᴴᵪᴳ′)
        _prune_at_threshold!(eachcol(tₐᵦ_basis_reim_ᴴᵪᴳ′_sparsified))

        # now we have the new terms - store them as `TightBindingTerm`s
        for tᴴᵪᴳ in eachcol(tₐᵦ_basis_reim_ᴴᵪᴳ′_sparsified)
            tbbᴴᵪᴳ = TightBindingBlock{D}(
                tbb.br1,
                tbb.br2,
                tbb.ordering1,
                tbb.ordering2,
                tbb.h_orbit,
                tbb.Mm,
                tᴴᵪᴳ,
                tbb.diagonal_block
            )
            h = TightBindingTerm{D}(
                tbt.axis,
                tbt.block_ij,
                tbbᴴᵪᴳ, #= .block =#
                tbt.hermiticity,
                tbt.brs,
            )
            push!(complement_tbs, h)
        end
    end
    return TightBindingModel{D}(complement_tbs, tbm.cbr, tbm.positions, tbm.N)
end

