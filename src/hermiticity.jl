
"""
    obtain_basis_free_parameters_hermiticity(
        h_orbit::HoppingOrbit{D},
        brₐ::BandRep{D},
        brᵦ::BandRep{D},
        antihermitian::Bool,
        orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
        orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
        Mm::AbstractMatrix{SparseMatrixCSC{Int, Int}} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ)
    ) where {D}

Constructs a basis for the coefficient vectors `t⁽ⁿ⁾` that span the space of Hermitian (or
`antihermitian` if `true`) TB Hamiltonians `Hₛₜ(k) = vᵢ(k) Mᵢⱼₛₜ tⱼ = vᵀ(k) M⁽ˢᵗ⁾ t`.
We do this by assuming that each coefficient vector `t` is sorted into a vector of the
form `[tᴿ; itᴵ]` so that we can take the complex conjugate by as `t* = σ₃t`, which can
then be moved onto `M⁽ˢᵗ⁾` instead of `t`. The constraint `Hₛₜ(k) = (H†)ₛₜ(k) = Hₜₛ*(k)`
can then be expressed as `vᵀ(k) M⁽ˢᵗ⁾ tⱼ = v*ᵀ(k) M⁽ᵗˢ⁾ σ₃ t = vᵀ(k) (Pᵀ M⁽ᵗˢ⁾σ₃) t`,
which requires that `t` be a solution to the nullspace `M⁽ˢᵗ⁾ - Pᵀ M⁽ᵗˢ⁾ σ₃ = 0`. We
cast this as `Z - Q = 0`, with `Z = M⁽ˢᵗ⁾` and `Q = Pᵀ M⁽ᵗˢ⁾ σ₃`.

## Notes
For anti-Hermitian symmetry, we require `Hₛₜ(k) = -Hₜₛ*(k)`, which translates to
`M⁽ˢᵗ⁾ + Pᵀ M⁽ᵗˢ⁾ σ₃ = 0`; i.e., simply swaps the sign of `Q`
"""
function obtain_basis_free_parameters_hermiticity(
    h_orbit::HoppingOrbit{D},
    brₐ::BandRep{D},
    brᵦ::BandRep{D},
    antihermitian::Bool,
    orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
    orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
    Mm::AbstractMatrix{SparseMatrixCSC{Int, Int}} = construct_M_matrix(h_orbit, brₐ, brᵦ, orderingₐ, orderingᵦ)
) where {D}
    # Determine a basis for coefficient vectors t⁽ⁿ⁾ that span the space of Hermitian (or
    # `antihermitian` if `true`) TB Hamiltonians Hₛₜ(k) = vᵢ(k) Mᵢⱼₛₜ tⱼ = vᵀ(k) M⁽ˢᵗ⁾ t. 
    # We do this by assuming that each coefficient vector `t` is sorted into a vector of the
    # form `[tᴿ; itᴵ]` so that we can take the complex conjugate by as `t* = σ₃t`, which can
    # then be moved onto `M⁽ˢᵗ⁾` instead of `t`. The constraint `Hₛₜ(k) = (H†)ₛₜ(k) = Hₜₛ*(k)` 
    # can then be expressed as `vᵀ(k) M⁽ˢᵗ⁾ tⱼ = v*ᵀ(k) M⁽ᵗˢ⁾ σ₃ t = vᵀ(k) (Pᵀ M⁽ᵗˢ⁾σ₃) t`,
    # which requires that `t` be a solution to the nullspace `M⁽ˢᵗ⁾ - Pᵀ M⁽ᵗˢ⁾ σ₃ = 0`. We
    # cast this as `Z - Q = 0`, with `Z = M⁽ˢᵗ⁾` and `Q = Pᵀ M⁽ᵗˢ⁾ σ₃`.
    # NB: For anti-Hermitian symmetry, we require `Hₛₜ(k) = -Hₜₛ*(k)`, which translates to
    #     `M⁽ˢᵗ⁾ + Pᵀ M⁽ᵗˢ⁾ σ₃ = 0`; i.e., simply swaps the sign of `Q`

    # Step 1: encode `M⁽ˢᵗ⁾` as the Z tensor (acting on `[tᴿ; itᴵ]`, i.e., as `[M⁽ˢᵗ⁾ M⁽ˢᵗ⁾]`)
    Z = [hcat(M, M) for M in Mm]

    # Step 2: encode Pᵀ M⁽ᵗˢ⁾σ₃ as the Q tensor
    Q = constraint_hermiticity(Mm, h_orbit, antihermitian)

    # Step 3: collate constraints across all `s` and `t`; matrix rows below are constraints
    constraints = _aggregate_constraints(Q, Z)

    # Step 4: find the nullspace of the constraints
    tₐᵦ_basis_matrix = nullspace(constraints; atol = NULLSPACE_ATOL_DEFAULT)

    return tₐᵦ_basis_matrix
end

function constraint_hermiticity(
    Mm::AbstractMatrix{<:SparseMatrixCSC{<:Number}},
    h_orbit::HoppingOrbit{D},
    antihermitian::Bool,
) where {D}
    opI = inversion(Val(D)) # inversion operation
    P = _permute_symmetry_related_hoppings_under_symmetry_operation(h_orbit, opI)
    Pᵀ = sparse(transpose(P))
    sgn = antihermitian ? -1 : 1 # anti-hermiticity flips signs
    return map(CartesianIndices(Mm)) do I
        s, t = Tuple(I)
        PᵀMᵗˢ = sgn * (Pᵀ * Mm[t, s]) # (±)Pᵀ M⁽ᵗˢ⁾
        hcat(PᵀMᵗˢ, -PᵀMᵗˢ)           # (±)Pᵀ M⁽ᵗˢ⁾ σ₃
    end
end
