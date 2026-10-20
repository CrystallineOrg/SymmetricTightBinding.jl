"""
    add_symmetry_constraints!(
        constraints::ConstraintRows,
        Mm::AbstractMatrix{SparseMatrixCSC{Int, Int}},
        h_orbit::HoppingOrbit{D},
        brₐ::BandRep{D},
        brᵦ::BandRep{D},
        gens::AbstractVector{<:AbstractOperation{D}},
       [orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
        orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ)]
    ) --> ConstraintRows

Add the constraints imposed on `t` by the spatial symmetries `gens` to `constraints`.

We have constructed the representation matrices such that gΦ(k) = ρᵀ(g)Φ(Rk); then, the
Hamiltonian transforms under symmetries as H(g𝐤) = ρₐₐ(g) H(𝐤) ρᵦᵦ⁺(g). The left-hand
side, H(g𝐤), is equivalent to permuting the elements of `v` (and hence the rows of each
`Mm[s,t]`; see `_permute_symmetry_related_hoppings_under_symmetry_operation`) as
Hₛₜ(g𝐤) = vᵀ(𝐤) Pᵀ M⁽ˢᵗ⁾ t. The right-hand side is

(ρₐₐ)ᵣₛ Hₛₜ (ρᵦᵦ⁻¹)ₜₗ = (ρₐₐ)ᵣₛ vᵢ Mᵢⱼₛₜ tⱼ (ρᵦᵦ⁻¹)ₜₗ = vᵢ (ρₐₐ)ᵣₛ Mᵢⱼₛₜ (ρᵦᵦ⁻¹)ₜₗ tⱼ,

i.e., vᵀ(𝐤) Q⁽ʳˡ⁾ t with Q⁽ʳˡ⁾ = Σₛₜ (ρₐₐ)ᵣₛ (ρᵦᵦ⁻¹)ₜₗ M⁽ˢᵗ⁾. Hence, `t` must satisfy
(Q⁽ʳˡ⁾ - Pᵀ M⁽ʳˡ⁾) t = 0 for every generator and every `(r,l)`; each row of
`Q⁽ʳˡ⁾ - Pᵀ M⁽ʳˡ⁾` is a complex constraint on `t`, added to `constraints` in realified form.

The sum over `s,t` runs only over the nonzero elements of `ρₐₐ` and `ρᵦᵦ`: as these are
block-permutation-like matrices, each Q⁽ʳˡ⁾ is a combination of only a few M⁽ˢᵗ⁾.
"""
function add_symmetry_constraints!(
    constraints::ConstraintRows,
    Mm::AbstractMatrix{SparseMatrixCSC{Int, Int}},
    h_orbit::HoppingOrbit{D},
    brₐ::BandRep{D},
    brᵦ::BandRep{D},
    gens::AbstractVector{<:AbstractOperation{D}},
    orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
    orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
) where {D}
    ρsₐₐ = site_induced_sgrep_excl_phase.(Ref(brₐ), gens, Ref(orderingₐ))
    ρsᵦᵦ = if brₐ == brᵦ # might as well reuse, if same BR
        ρsₐₐ
    else
        site_induced_sgrep_excl_phase.(Ref(brᵦ), gens, Ref(orderingᵦ))
    end

    Q_terms = ComplexConstraintTerm{SparseMatrixCSC{Int, Int}}[] # reused for each block
    for (op, ρₐₐ, ρᵦᵦ) in zip(gens, ρsₐₐ, ρsᵦᵦ)
        Pᵀ = _permute_symmetry_related_hoppings_under_symmetry_operation(h_orbit, op)
        for l in axes(Mm, 2), r in axes(Mm, 1)
            # Q⁽ʳˡ⁾ - Pᵀ M⁽ʳˡ⁾, with Q⁽ʳˡ⁾ = Σₛₜ (ρₐₐ)ᵣₛ (ρᵦᵦ)*ₗₜ M⁽ˢᵗ⁾ (nonzero terms only)
            empty!(Q_terms)
            append!(Q_terms, (ComplexConstraintTerm(ρₐₐ[r, s] * conj(ρᵦᵦ[l, t]), Mm[s, t])
                              for s in axes(ρₐₐ, 2) if !iszero(ρₐₐ[r, s])
                              for t in axes(ρᵦᵦ, 2) if !iszero(ρᵦᵦ[l, t])))
            add_complex_constraints!(constraints, Q_terms,
                                     (ComplexConstraintTerm(-1, Pᵀ * Mm[r, l]),))
        end
    end
    return constraints
end

"""
Build the P matrix for a particular symmetry operation acting on k-space, whose transpose
`Pᵀ` permutes the rows of the M matrix; `Pᵀ` (not `P`) is returned, as a `RowPermutation`.

To obtain the P matrix, we exploit that the action is on exponentials of the type
``exp(-2πi𝐤⋅δ)``, and instead act on δ ∈ `h_orbit.orbit` rather than on k. Because of this,
we need to use the inverse of the rotation part of the symmetry operation.

!!! details "Sketch of proof"
    Assume g={R|τ} and `Crystalline` implements gk=(R⁻¹)ᵀk. Then
    (gk)⋅δ = ((R⁻¹)ᵀk)⋅δ + τ = k⋅(R⁻¹)δ.

!!! info
    It is assumed that the operation `op` is provided in a primitive setting.

!!! warning
    This function is an internal helper function for `add_symmetry_constraints!` and
    `add_hermiticity_constraints!` and is not part of the public API.
"""
function _permute_symmetry_related_hoppings_under_symmetry_operation(
    h_orbit::HoppingOrbit{D},
    op::AbstractOperation{D}, # only the rotation part is used
) where {D}
    # P is a square matrix that acts as `op` on `v`, i.e., represents ``op ∘ v``, via a
    # matrix-vector product `P*v`. Equivalently, `P` can act on `M` (via its transpose)
    # for expressions of the kind (Pv)ᵀ M … = vᵀ Pᵀ M …
    
    N = length(orbit(h_orbit))
    cols = Vector{Int}(undef, N) # `P[i, cols[i]] = 1`, i.e., the nonzero column of each row
    for (i, δᵢ) in enumerate(orbit(h_orbit))
        # crystalline implements gk = [R⁻¹]ᵀk. However, since we don't have access to k,
        # being a symbolic variable, we will translate its action into δ. For this,
        # we implement a simple trick: ([R⁻¹]ᵀ k)·r = R⁻¹ᵢⱼ kᵢ rⱼ = k·(R⁻¹ r)
        R⁻¹ = SymOperation(inv(rotation(op))) # rotation-only part of `op`
        δᵢ′ = compose(R⁻¹, δᵢ)
        j = findfirst(orbit(h_orbit)) do δ′′
            isapprox(δᵢ′, δ′′, nothing, false; atol = VEC_CMP_ATOL)
        end
        isnothing(j) &&
            error(lazy"hopping element $δᵢ not closed under $op in $(orbit(h_orbit))")
        cols[i] = j
        # since P acts as g v = P v, we have (gv)ᵢ = vⱼ = ∑ₖ Pᵢₖ vₖ so Pᵢₖ = δᵢⱼ
    end
    # since `Π[σ[j], j] = 1` for a `Π::RowPermutation`, `cols` directly gives `Pᵀ`, which
    # is what callers need (`Pᵀ[cols[i], i] = P[i, cols[i]] = 1`)
    Pᵀ = RowPermutation(cols)
    return Pᵀ
end
