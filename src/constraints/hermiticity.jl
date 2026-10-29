"""
    add_hermiticity_constraints!(
        constraints::ConstraintRows,
        Mm::AbstractMatrix{SparseMatrixCSC{Int, Int}},
        h_orbit::HoppingOrbit{D},
        antihermitian::Bool,
    ) --> ConstraintRows

Add the constraints on the coefficient vectors `t` that make a (diagonal) block of the TB
Hamiltonian `Hₛₜ(k) = vᵢ(k) Mᵢⱼₛₜ tⱼ = vᵀ(k) M⁽ˢᵗ⁾ t` Hermitian (or `antihermitian` if
`true`) to `constraints`.
We do this by assuming that each coefficient vector `t` is sorted into a vector of the
form `[tᴿ; itᴵ]` so that we can take the complex conjugate by as `t* = σ₃t`, which can
then be moved onto `M⁽ˢᵗ⁾` instead of `t`. The constraint `Hₛₜ(k) = (H†)ₛₜ(k) = Hₜₛ*(k)`
can then be expressed as `vᵀ(k) M⁽ˢᵗ⁾ tⱼ = v*ᵀ(k) M⁽ᵗˢ⁾ σ₃ t = vᵀ(k) (Pᵀ M⁽ᵗˢ⁾σ₃) t`,
which requires that `t` be a solution to the nullspace `Pᵀ M⁽ᵗˢ⁾ σ₃ - M⁽ˢᵗ⁾ = 0`, with
`M⁽ˢᵗ⁾` acting on `[tᴿ; tᴵ]` as `[M⁽ˢᵗ⁾ M⁽ˢᵗ⁾]`.

## Notes
For anti-Hermitian symmetry, we require `Hₛₜ(k) = -Hₜₛ*(k)`, which translates to
`Pᵀ M⁽ᵗˢ⁾ σ₃ + M⁽ˢᵗ⁾ = 0`; i.e., simply swaps the sign of the `Pᵀ M⁽ᵗˢ⁾ σ₃` term.
"""
function add_hermiticity_constraints!(
    constraints::ConstraintRows,
    Mm::AbstractMatrix{SparseMatrixCSC{Int, Int}},
    h_orbit::HoppingOrbit{D},
    antihermitian::Bool,
) where {D}
    # `P` represents the action of inversion on `v`, i.e., `v*(k) = v(-k) = P v(k)`
    opI = inversion(Val(D))
    Pᵀ = _permute_symmetry_related_hoppings_under_symmetry_operation(h_orbit, opI)
    sgn = antihermitian ? -1 : 1 # anti-hermiticity flips signs

    for t in axes(Mm, 2), s in axes(Mm, 1)
        # (±)Pᵀ M⁽ᵗˢ⁾ σ₃ - [M⁽ˢᵗ⁾ M⁽ˢᵗ⁾], acting on `[tᴿ; tᴵ]`
        add_constraints!(constraints, (RealConstraintTerm(-1, -1, Mm[s, t]),),
                                      (RealConstraintTerm(sgn, -sgn, Pᵀ * Mm[t, s]),))
    end
    return constraints
end
