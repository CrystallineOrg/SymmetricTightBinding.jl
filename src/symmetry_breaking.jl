"""
    subduced_complement(tbm::TightBindingModel{D}, Rs, sgnumᴴ::Int; timereversal)
                                                        --> TightBindingModel{D}

Given a model `tbm` associated with a space group ``G``, determine the new, independent
tight-binding terms (i.e., the orthogonal complement of terms) that become
symmetry-allowed when the model's space group is reduced to a subgroup ``H ≤ G`` with space
group number `sgnumᴴ` and time-reversal symmetry `timereversal`.

Practically, the function answers the question: which new tight-binding terms become allowed
if the symmetry of the model is reduced from space group ``G`` to subgroup ``H``?

`Rs` is the hopping range that `tbm` was built with, i.e., `tbm = tb_hamiltonian(cbr, Rs)`.

!!! note "Why is `Rs` needed?"
    `tbm` may lack terms for some hoppings in `Rs`, namely those forbidden in ``G``; these
    may nevertheless be allowed in ``H``.

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

julia> brs = bandreps(17, Val(2); timereversal = true);

julia> cbr = @composite brs[5]

julia> Rs = [[0,0], [1,0]];

julia> tbm = tb_hamiltonian(cbr, Rs)
```
Each of the 4 terms in this model is proportional to an identity matrix at K = (1/3, 1/3).
Using `subduced_complement`, we can find the new terms that appear if we imagine lowering
the symmetry from plane group ⋕17 to ⋕16 (which has no mirror symmetry) while also removing
time-reversal symmetry.
```julia-repl
julia> Δtbm = subduced_complement(tbm, Rs, 16; timereversal = false)
2-term 2×2 TightBindingModel{2} (hermitian, spinless) over (2b|A₁), where zᵢ=exp(-2πik·δᵢ):
┌─
1. ⎡ -iz₁+iz̄₁-iz₂+iz̄₂+iz₃-iz̄₃  0                       ⎤
│  ⎣ 0                         iz₁-iz̄₁+iz₂-iz̄₂-iz₃+iz̄₃ ⎦
└─ (2b|A₁) self-term.  δ₁=[1,0], δ₂=[0,1], δ₃=[1,1]
┌─
2. ⎡ 0                  z̄₁+z̄₂+z₃-z̄₄-z₅-z̄₆ ⎤
│  ⎣ z₁+z₂+z̄₃-z₄-z̄₅-z₆  0                 ⎦
└─ (2b|A₁) self-term.  δ₁=[4/3,-1/3], δ₂=[1/3,5/3], δ₃=[5/3,4/3], δ₄=[1/3,-4/3], δ₅=[5/3,1/3], δ₆=[4/3,5/3]
```
The first of the of these terms is not diagonal at K and so opens a gap at the Dirac point:
```julia-repl
julia> Δtbm[1](ReciprocalPoint(1/3, 1/3))
2×2 Matrix{ComplexF64}:
 -5.19615+2.77556e-16im      0.0+0.0im
      0.0+0.0im          5.19615+4.996e-16im
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
function subduced_complement(
    tbm::TightBindingModel{D},
    Rs::AbstractVector{<:AbstractVector{<:Integer}},
    sgnumᴴ::Int;
    kws...
) where D
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
    pᴴ²ᴳ = -Pᴴ²ᴳ * t.p

    _gensᴴ = generators(sgnumᴴ, SpaceGroup{D}) # in H setting
    gensᴴ = transform.(_gensᴴ, Ref(Pᴴ²ᴳ), Ref(pᴴ²ᴳ))
    if isspinful(tbm)
        # attach the SU(2) elements of G's setting; transforming H's double group generators
        # instead would keep the SU(2) elements of H's setting (cf. Crystalline's `SU2`)
        # (Ē, the only extra generator of H's double group, imposes no constraint)
        gensᴴ = [DSymOperation{D}(g, SU2(g, sgnumᴳ)) for g in gensᴴ]
    end

    return _subduced_complement(tbm, Rs, gensᴴ; kws...)
end

"""
    _subduced_complement(
        tbm::TightBindingModel{D},
        Rs,
        gensᴴ::AbstractVector{<:AbstractOperation{D}};
        timereversal::Bool
    ) where D --> TightBindingModel{D}

Implementation of [`subduced_complement`](@ref), taking the generators `gensᴴ` of the
subgroup ``H`` rather than its space group number.

`gensᴴ` must be given in the *conventional* setting of the original group ``G`` (i.e., in
the setting of `tbm`); they are converted to the primitive setting internally. This is why
the method is not part of the public API: obtaining `gensᴴ` in G's setting requires the
transformation dance of the `sgnumᴴ` method, which is not reasonable to ask of a caller.

For spinful models, ``G`` and ``H`` refer to the double group, and `gensᴴ` must accordingly
be double group operations (`DSymOperation`).

!!! warning
    This function is an internal helper function for `subduced_complement` and is not part
    of the public API.
"""
function _subduced_complement(
    tbm::TightBindingModel{D, S, IR, SIR},
    Rs::AbstractVector{<:AbstractVector{<:Integer}},
    gensᴴ::AbstractVector{<:AbstractOperation{D}};
    timereversal::Bool = first(tbm.cbr.brs).timereversal, # ← whether H has time-reversal
) where {D, S, IR, SIR}
    _check_operation_spin(tbm, gensᴴ) # check `tbm` & `gensᴴ` have equal `isspinful`
    timereversalᴳ = first(tbm.cbr.brs).timereversal
    if timereversalᴳ == false && timereversal == true
        error(
            "requested subgroup `timereversal = true`, but original model was built without time-reversal present; input for H must maintain or reduce symmetry",
        )
    end

    sgnumᴳ = num(tbm.cbr)
    @assert _issubgroup(gensᴴ, sgnumᴳ) LazyString(
        "`gensᴴ` must generate a subgroup of G (⋕", sgnumᴳ, ") and be given in G's \
         conventional setting, but ", gensᴴ, " is not a subset of the operations of G")

    # the constraint machinery in `_obtain_basis_free_parameters` works in the primitive
    # setting (cf. `obtain_basis_free_parameters`), but `gensᴴ` is given in the conventional
    # setting of G: convert, lest we compare conventional-setting operations against the
    # primitivized site symmetry groups of `site_induced_sgrep` (which finds no
    # matching coset and errors out)
    gensᴴ′ = primitivized_generators(gensᴴ, sgnumᴳ)

    # each `tbm[i]` term lives in one block of the Hamiltonian and on one hopping orbit; 
    # each such (block, orbit) pair has its own coefficient basis, so we compute the
    # complement pair by pair: `groups` holds a `(; block_ij, block, idxs)` for each pair
    # over `Rs`, with `block` a representative block and `idxs` the indices of the terms of
    # `tbm` on the pair (empty if G allows no term on it)
    groups = _subduction_groups(tbm, Rs)
    complement_tbs = TightBindingTerm{D, S, IR, SIR}[]
    isempty(groups) && return TightBindingModel(complement_tbs, tbm.cbr, tbm.positions, tbm.N)
    axis, brs = first(tbm.terms).axis, first(tbm.terms).brs # shared by all terms of `tbm`

    # now we can compute a new coefficient basis in H and compare with our original basis,
    # progressing "group by group"
    for (; block_ij, block, idxs) in groups
        # first, compute basis of coefficients for new subset of generators (`gensᴴ`)
        tₐᵦ_basis_reimᴴ_vs = _obtain_basis_free_parameters(
            block.h_orbit,
            block.br1,
            block.br2,
            block.ordering1,
            block.ordering2,
            block.Mm,
            gensᴴ′,
            timereversal,
            block_ij[1] == block_ij[2], #= .diagonal_block =#
            S,                          #= hermiticity =#
        )
        # check output dimensions
        if length(tₐᵦ_basis_reimᴴ_vs) < length(idxs)
            error(
                LazyString(
                    "unexpectedly found lower-dimensional basis space for model in \
                     subduced group (dim ",
                    length(tₐᵦ_basis_reimᴴ_vs),
                    " < ",
                    length(idxs),
                    " terms, for block ",
                    block_ij,
                    " & orbit ",
                    representative(block.h_orbit),
                    "); unexpected and unhandled - make sure the generators are a subset \
                     of the original generators (i.e., that fewer constraints apply than \
                     originally)",
                ),
            )
        elseif length(tₐᵦ_basis_reimᴴ_vs) == length(idxs)
            continue # basis must then be unchanged; nothing to add for this index group
        end

        if isempty(idxs)
            # nothing is spanned in G, so the entire H basis is the complement
            # it is already in the same sparsified form that the projection & SVD below would
            # return it in, so we can store the terms directly
            for tᴴ in tₐᵦ_basis_reimᴴ_vs
                tbbᴴ = TightBindingBlock{D, S}(
                    block.br1,
                    block.br2,
                    block.ordering1,
                    block.ordering2,
                    block.h_orbit,
                    block.Mm,
                    tᴴ,
                    block.diagonal_block
                )
                h = TightBindingTerm(axis, block_ij, tbbᴴ, brs)
                push!(complement_tbs, h)
            end
            continue
        end

        # get "original" coefficient basis in G from `tbm[idxs]`
        tₐᵦ_basis_reimᴴ = stack(tₐᵦ_basis_reimᴴ_vs)
        tₐᵦ_basis_reimᴳ = Matrix{Float64}(undef, length(block.t), length(idxs))
        for (n, i) in enumerate(idxs)
            tbbᵢ = tbm.terms[i].block
            tₐᵦ_basis_reimᴳ[:, n] .= tbbᵢ.t
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
        tₐᵦ_basis_reim_ᴴᵪᴳ′ = Matrix{Float64}(undef, length(block.t), Nᴴ)
        for (n, (u, σ)) in enumerate(zip(eachcol(Uᴴᵪᴳ), σs))
            n > Nᴴ && continue # there should be exactly Nᴴ non-zero singular values
            if σ < NULLSPACE_ATOL_DEFAULT
                # make sure we don't have any surprises & verify our assumption above
                error(
                    LazyString(
                        "unexpectedly found near-zero singular value for a SVD \
             column vector that ought to have been a proper basis vector (Nᴴ = ",
                        Nᴴ,
                        " σs = ",
                        σs,
                        ", for block ",
                        block_ij,
                        " & orbit ",
                        representative(block.h_orbit),
                        ")",
                    ),
                )
            end
            # keep this vector
            tₐᵦ_basis_reim_ᴴᵪᴳ′[:, n] .= u
        end

        # finally, make the basis vectors we have now "pretty"
        tₐᵦ_basis_reim_ᴴᵪᴳ′_sparsified = _poormans_sparsification(tₐᵦ_basis_reim_ᴴᵪᴳ′)
        _prune_at_threshold!(eachcol(tₐᵦ_basis_reim_ᴴᵪᴳ′_sparsified))

        # now we have the new terms - store them as `TightBindingTerm`s
        for tᴴᵪᴳ in eachcol(tₐᵦ_basis_reim_ᴴᵪᴳ′_sparsified)
            tbbᴴᵪᴳ = TightBindingBlock{D, S}(
                block.br1,
                block.br2,
                block.ordering1,
                block.ordering2,
                block.h_orbit,
                block.Mm,
                tᴴᵪᴳ,
                block.diagonal_block
            )
            h = TightBindingTerm(axis, block_ij, tbbᴴᵪᴳ, brs)
            push!(complement_tbs, h)
        end
    end
    return TightBindingModel(complement_tbs, tbm.cbr, tbm.positions, tbm.N)
end

"""
    _issubgroup(gensᴴ::AbstractVector{<:AbstractOperation{D}}, sgnumᴳ::Int)  --> Bool

Return whether every operation of `gensᴴ` is an operation of the space group `sgnumᴳ` (up to
lattice translations), i.e., whether `gensᴴ` generates a subgroup ``H ≤ G``.

`gensᴴ` is assumed given in the conventional setting of ``G``. For double group operations
(`DSymOperation`s), the comparison is with the double group of ``G``, including its SU(2)
elements.

!!! warning
    This function is an internal helper function for `subduced_complement` and is not part
    of the public API.
"""
function _issubgroup(
    gensᴴ::AbstractVector{O}, sgnumᴳ::Int
) where {D, O<:AbstractOperation{D}}
    cntr = centering(sgnumᴳ, D)
    opsᴳ = spacegroup(sgnumᴳ, Val(D); spinful = Val(isspinful(O)))
    return all(opᴴ -> any(opᴳ -> isapprox(opᴳ, opᴴ, cntr), opsᴳ), gensᴴ)
end

"""
    _subduction_groups(tbm::TightBindingModel{D, S}, Rs)
                        --> Vector{@NamedTuple{block_ij, block, idxs}}

The (block, orbit) pairs that `subduced_complement` must visit, i.e., those generated
internally in `tb_hamiltonian(tbm.cbr, Rs)`. Each is given by its block index `block_ij`, a
representative `block::TightBindingBlock` - from which the orbit and M-tensor are read - and
the indices `idxs` of the terms of `tbm` that span its coefficient basis in ``G``.

A pair may carry no term in `tbm` (empty `idxs`): if every coefficient is forbidden in ``G``,
the pair is then represented by a zero-coefficient block (issue #117). Terms of `tbm` on
orbits not generated by `Rs` are ignored.

!!! warning
    This function is an internal helper function for `subduced_complement` and is not part
    of the public API.
"""
function _subduction_groups(
    tbm::TightBindingModel{D, S, IR, SIR},
    Rs
) where {D, S, IR, SIR}
    groups = @NamedTuple{block_ij::NTuple{2, Int},
                         block::TightBindingBlock{D, S, IR, SIR},
                         idxs::Vector{Int}}[]
    isempty(tbm) && return groups

    # the block structure of the model, exactly as `tb_hamiltonian` built it
    for block_info in _block_iterates(first(tbm).brs, Rs, Val(S))
        (; block_ij, br1, br2, ordering1, ordering2, diagonal_block, h_orbits) = block_info
        for h_orbit in h_orbits
            idxs = findall(tbm) do tbt
                tbt.block_ij == block_ij &&
                    isapproxin(representative(h_orbit), orbit(tbt.block.h_orbit),
                               nothing, false; atol = VEC_CMP_ATOL)
            end
            block = if !isempty(idxs)
                # use a block of `tbm` itself: its coefficients are indexed by its own orbit,
                # which a re-enumeration only reproduces as a set, not in order
                tbm.terms[first(idxs)].block
            else
                # a pair with no G-allowed term: stand it in with a `t = zeros(…)` block
                Mm = construct_M_matrix(h_orbit, br1, br2, ordering1, ordering2)
                TightBindingBlock{D, S}(br1, br2, ordering1, ordering2, h_orbit, Mm,
                                        #=t=# zeros(2size(Mm, 2)), diagonal_block)
            end
            push!(groups, (; block_ij, block, idxs))
        end
    end
    return groups
end
