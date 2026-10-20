#=
** Implementation and theory notes **

TRS can be understood as a spatial symmetry when acting on the Hamiltonian:
   D(𝒯)H(k)D(𝒯)† = H(𝒯k) → Γ(𝒯)H*(k)Γ(𝒯)† = H(-k),
where D is the whole operator and Γ is only its unitary part, so D(𝒯) = Γ(𝒯)𝒯.
For *spinless* band representations, 𝒯² = +1, the basis of the representation is real and
Γ(𝒯) can be chosen as Γ(𝒯) = I so that H*(k) = H(-k).
For *spinful* band representations, 𝒯² = -1 and we can (and do, cf. Crystalline's
`timereversal_unitary`) choose Γ(𝒯) = J ⊗ 𝟙 with J = iσʸ ⊗ 𝟙ₙ (J over the partner
functions at each site: σʸ over the Kramers index, 𝟙ₙ over the remaining n of a 2n-dim.
site-symmetry irrep; and 𝟙 over the sites of the Wyckoff position's orbit - overall,
following the partner-function-major order of `OrbitalOrdering`). In both cases, Γ is
real (a signed permutation).

NB: The "realification" of the representation, i.e., the basis choice that brings Γ to
the above forms, is performed in Crystalline via `physical_realify`, applied when
`timereversal = true`. This is done in `bandreps` itself, i.e., can be assumed here.

We have that Hₛₜ(k) = vᵀ(k) Mₛₜ t, with s and t indexing the orbitals of `brₐ` and `brᵦ`,
and Mₛₜ = `Mm[s,t]` (cf. `construct_M_matrix`). Remember that we split t into real and
imaginary parts: t = [real(t); i imag(t)], and Mm acts as [Mm Mm]. Thus, we can rewrite the
Hamiltonian as:
   Hₛₜ(k) = vᵀ(k) [Mₛₜ Mₛₜ] [real(t); i imag(t)]

Thus, the action of TRS on the Hamiltonian, H(-k) = ΓH*(k)Γ†, can be decomposed into parts:

1. Hₛₜ(-k) = vᵀ(-k) [Mₛₜ Mₛₜ] [real(t); i*imag(t)]

2. (ΓH*(k)Γ†)ₛₜ = Σₘₙ Γₐ,ₛₘ H*ₘₙ(k) Γᵦ,ₜₙ                  (NB: Γ is real, so Γ† = Γᵀ)
                = Σₘₙ Γₐ,ₛₘ (v*)ᵀ(k) [Mₘₙ Mₘₙ] [real(t); -i*imag(t)] Γᵦ,ₜₙ
                = (v*)ᵀ(k) (Σₘₙ Γₐ,ₛₘ [Mₘₙ Mₘₙ] Γᵦ,ₜₙ) [real(t); -i*imag(t)]
                = (v*)ᵀ(k) [M̃ₛₜ -M̃ₛₜ] [real(t); i*imag(t)]
                = vᵀ(-k) [M̃ₛₜ -M̃ₛₜ] [real(t); i*imag(t)]
   with M̃ₛₜ = Σₘₙ Γₐ,ₛₘ Mₘₙ Γᵦ,ₜₙ, using that Γ is real and that v*(k) = v(-k).
   In practice, we compute M̃ as the product M̃ = Γₐ Mm Γᵦᵀ over the orbital (s,t) axes of
   Mm, i.e., as M̃ₛₜ = Σₘₙ Γₐ,ₛₘ Γᵦ,ₜₙ Mm[m,n].

Imposing the condition H(-k) = ΓH*(k)Γ† we get:
  vᵀ(-k) [Mₛₜ-M̃ₛₜ Mₛₜ+M̃ₛₜ] [real(t); i*imag(t)] = 0   (for all s, t)
For spinless band representations, M̃ = M, and this is simply [0 2Mₛₜ] (we use
[0 Mₛₜ], which has the same null space). For spinful band representations, M̃ permutes
the orbitals of each site among Kramers partners (with signs), so TRS relates hoppings,
e.g., t↑↑ = t↓↓* and t↑↓ = -t↓↑*. Note that Γ enters twice, so its overall sign is
irrelevant.

This way of casting the problem is very convenient (and doesn't require any modification of 
the v vectors to include "reversed" hoppings, which is not in general a necessary feature).
We just need to add the rows of [Mₛₜ-M̃ₛₜ Mₛₜ+M̃ₛₜ] for all s, t to the constraints on the
re/im-doubled coefficient vector (cf. `add_constraints!`). NB: each such row, taken on its
own, mixes the conditions on real(t) and imag(t); but the rows for (s,t) and its Kramers
partner (s̄,t̄) (with M̃ₛ̄ₜ̄ = ±Mₛₜ) combine to separate conditions on real(t) and imag(t).
=#

"""
    add_timereversal_constraints!(
        constraints::ConstraintRows,
        Mm::AbstractMatrix{SparseMatrixCSC{Int, Int}},
        brₐ::BandRep{D},
        brᵦ::BandRep{D},
       [orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
        orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ)]
    ) --> ConstraintRows

Add the constraints imposed on `t` by time-reversal symmetry to `constraints`: i.e., the
rows of [Mₛₜ-M̃ₛₜ Mₛₜ+M̃ₛₜ] for all s, t, acting on `[tᴿ; tᴵ]` (see notes above).
"""
function add_timereversal_constraints!(
    constraints::ConstraintRows,
    Mm::AbstractMatrix{SparseMatrixCSC{Int, Int}},
    brₐ::BandRep{D},
    brᵦ::BandRep{D},
    orderingₐ::OrbitalOrdering{D} = OrbitalOrdering(brₐ),
    orderingᵦ::OrbitalOrdering{D} = OrbitalOrdering(brᵦ),
) where {D}
    S = isspinful(brₐ)
    S == isspinful(brᵦ) || error("input `BandRep`s must have the same spin")

    if !S # spinless: `M̃ = M`, so the constraint rows are `[0 M]`
        for M in Mm
            add_constraints!(constraints, (RealConstraintTerm(0, 1, M),))
        end
    else
        # spinful: the constraint rows are `[M-M̃ M+M̃]`, with `M̃ = Γₐ M Γᵦᵀ` acting on
        # the orbital axes (`Γₐ` & `Γᵦ`: unitary parts of time reversal across the orbitals
        # of `brₐ` & `brᵦ`, in the orbital ordering of `OrbitalOrdering`)
        Γₐ = site_induced_timereversal_unitary(brₐ, orderingₐ) # defined here only, since they
        Γᵦ = site_induced_timereversal_unitary(brᵦ, orderingᵦ) # aren't needed (≡1) if spinless
        terms = RealConstraintTerm{SparseMatrixCSC{Int, Int}}[] # reused for each block
        for t in axes(Mm, 2), s in axes(Mm, 1)
            # [M⁽ˢᵗ⁾-M̃⁽ˢᵗ⁾  M⁽ˢᵗ⁾+M̃⁽ˢᵗ⁾], with M̃⁽ˢᵗ⁾ = Σₘₙ (Γₐ)ₛₘ (Γᵦ)ₜₙ M⁽ᵐⁿ⁾; since `Γₐ` &
            # `Γᵦ` are signed permutations, M̃⁽ˢᵗ⁾ involves only a few M⁽ᵐⁿ⁾
            empty!(terms)
            push!(terms, RealConstraintTerm(1, 1, Mm[s, t]))
            for m in axes(Γₐ, 2), n in axes(Γᵦ, 2)
                γ = Γₐ[s, m] * Γᵦ[t, n]
                iszero(γ) || push!(terms, RealConstraintTerm(-γ, γ, Mm[m, n]))
            end
            add_constraints!(constraints, terms)
        end
    end
    return constraints
end
