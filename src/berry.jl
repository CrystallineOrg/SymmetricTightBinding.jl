"""
    berrycurvature(
        ptbm::ParameterizedTightBindingModel{D},
        k::ReciprocalPointLike{D},
        n::Integer,
        [∇Hs::NTuple{D, Matrix{ComplexF64}}]
    )  --> Ω

Compute the Berry curvature `Ω` of the `n`th band of the coefficient-parameterized
tight-binding model `ptbm` at the **k**-point `k`.
If `D == 3`, `Ω` is a 3-element vector, representing, the three components of the Berry
curvature pseudovector ``[Ω₁, Ω₂, Ω₃] = [Ω₂₃, Ω₃₁, Ω₁₂]``. 
If `D == 2`, `Ω` is a scalar, representing the Berry curvature pseudoscalar ``Ω₃ = Ω₁₂``.

The optional `∇Hs` argument is a tuple of work-matrices used to store the momentum gradient
of the Hamiltonian.

## Implementation

The Berry curvature is evaluated using the Kubo-like or "sum-over-states" formula [^1]:
```math
\\Omega_{\\mu, \\nu}^n
=
\\sum_{n' \\neq n}
\\frac{\\langle n |\\partial H/\\partial k_\\mu|n'\\rangle \\langle n'|\\partial H/\\partial k_\\nu|n>
       - (\\nu \\leftrightarrow \\mu)}
      {(E_n - E_{n'})^2}
=
-2\\mathrm{Im} \\sum_{n' \\neq n}
\\frac{\\langle n |\\partial H/\\partial k_\\mu|n'\\rangle \\langle n'|\\partial H/\\partial k_\\nu|n>}
      {(E_n - E_{n'})^2}
```
with ``\\Omega_{\\mu, \\nu}^n`` denoting  the Berry curvature tensor component of the `n`th
band and ``E_n`` denoting the energy of the `n`th band. The Berry curvature pseudovector and
(antisymmetric) tensor components are related by ``\\Omega_1 = \\Omega_{23}``,
``\\Omega_2 = \\Omega_{31}``, and ``\\Omega_3 = \\Omega_{12}``.

[1]: Xiao, Chang, & Niu, "Berry phase effects on electronic properties", 
    [Rev. Mod. Phys. **82**, 1959 (2010)](https://doi.org/10.1103/RevModPhys.82.1959) (see
    Eq. (1.13)).

## Extended help

As elsewhere in SymmetricTightBinding.jl, this function operates in *reduced* coordinates,
with the momentum `k = (k₁, k₂, k₃)` indicating coefficients relative to the primitive
reciprocal lattice vectors ``\\mathbf{b}_i``, i.e., ``\\mathbf{k} = k_1\\mathbf{b}_1 +
k_2\\mathbf{b}_2 + \\ldots``.
Crucially, this also implies that returned Berry curvature components are *not* the
Cartesian components, but components relative to the primitive reciprocal basis.

**In 3D**, this implies that the returned vector `Ω = [Ω₁, Ω₂, Ω₃]` relates to the Berry
curvature pseudovector ``\\boldsymbol{\\Omega}`` via
```math
\\boldsymbol{\\Omega}
=
\\Omega_1 \\mathbf{b}_1 + \\Omega_2 \\mathbf{b}_2 + \\Omega_3 \\mathbf{b}_3
=
\\Omega_{23} \\mathbf{b}_1 + \\Omega_{31} \\mathbf{b}_2 + \\Omega_{12} \\mathbf{b}_3
```
Here, the tensor components ``\\Omega_{12}`` etc. correspond to momentum derivatives in the
*reduced* coordinates, i.e., the derivatives ``\\partial/\\partial k_{1,2}`` are taken with
respect to  reduced (dimensionless) coordinates ``k_{1,2}`` rather than Cartesian ones.

As a result, the returned Berry curvature pseudovector components ``[Ω₁, Ω₂, Ω₃]`` are
related to the Cartesian pseudovector components ``[Ω_x, Ω_y, Ω_z]`` by the primitive
reciprocal lattice matrix ``\\mathbf{B} = [\\mathbf{b}_1 \\mathbf{b}_2 \\mathbf{b}_3]``:
```math
[Ω_x, Ω_y, Ω_z]^{\\text{T}} = \\mathbf{B} [Ω₁, Ω₂, Ω₃]^{\\text{T}}
```

**In 2D**, the returned Berry curvature is a (pseudo)scalar ``Ω₃ = Ω₁₂``. This value
represents the component of the Berry curvature 2-form in the reduced coordinate system. It
is related to the physical pseudoscalar ``\\Omega_{xy}`` by the signed area of the Brillouin
zone ``A_{\\text{BZ}} = (\\mathbf{b}_1 \\times \\mathbf{b}_2) \\cdot \\hat{\\mathbf{z}}``:
```math
\\Omega_z = \\Omega_{xy} = \\Omega_{12} / A_{\\text{BZ}}
```
While this distinction is important for interpreting the curvature at an individual
**k**-point, it can be ignored when calculating the Chern number as an integral over the
reduced coordinates ``k_{1,2} \\in (-1/2, 1/2]``, since the geometric factors from the
curvature and the area element ``\\mathrm{d}^2\\mathbf{k}`` cancel out.

### Surface integrals and flux in 3D

When using the returned 3D Berry curvature pseudovector `[Ω₁, Ω₂, Ω₃]` to calculate
physical observables like the Chern number of a slice or the flux through a surface, it's
crucial to integrate the correct physical quantity. The assumed goal is always to compute
a flux integral:
```math
\\Phi = \\int_S \\boldsymbol{\\Omega} \\cdot \\,\\mathrm{d}\\mathbf{S}.
````

While this function returns the components `[Ω₁, Ω₂, Ω₃]` in the reduced basis, the flux
integrand requires the full dot product between the physical vector ``\\boldsymbol{\\Omega}
= \\Omega_1 \\mathbf{b}_1 + \\Omega_2 \\mathbf{b}_2 + \\Omega_3 \\mathbf{b}_3`` and the
physical area element ``\\mathrm{d}\\mathbf{S}``.

The implications of this statement can be appreciated if we consider e.g., calculating the
Chern number on a planar slice of the BZ, for example the plane spanned by ``\\mathbf{b}_1``
and ``\\mathbf{b}_2``. The physical area element is ``\\mathrm{d}\\mathbf{S} = 
(\\mathbf{b}_1 \\times \\mathbf{b}_2) \\, \\mathrm{d}k_1 \\mathrm{d}k_2``, so that the flux
integrand is:

```math
\\mathrm{d}\\Phi
=
\\boldsymbol{\\Omega} \\cdot d\\mathbf{S}
=
(\\Omega_1 \\mathbf{b}_1 + \\Omega_2 \\mathbf{b}_2 + \\Omega_3 \\mathbf{b}_3)
\\cdot
(\\mathbf{b}_1 \\times \\mathbf{b}_2) \\, \\mathrm{d}k_1 \\mathrm{d}k_2
```

By the properties of the triple product, the terms proportional to ``\\Omega_{1,2}`` vanish,
leaving:

```math
\\mathrm{d}\\Phi
=
\\Omega_3 (\\mathbf{b}_3 \\cdot (\\mathbf{b}_1 \\times \\mathbf{b}_2))
\\, \\mathrm{d}k_1 \\mathrm{d}k_2
```

The term `\\mathbf{b}_3 \\cdot (\\mathbf{b}_1 \\times \\mathbf{b}_2)` is the volume of the
3D Brillouin zone, `V_{\text{BZ}}`.

This shows that the physical flux through the ``(k₁, k₂)`` plane, involves the returned 3D
component `\\Omega_3``, but scaled by the 3D cell volume `V_BZ`.

For any general surface, the same principle applies: the full dot product must be evaluated,
which automatically incorporates the required geometric scaling factors.
"""
function berrycurvature(
    ptbm::ParameterizedTightBindingModel{D},
    k::ReciprocalPointLike{D},
    n::Integer,
    ∇Hs::NTuple{D, Matrix{ComplexF64}} = begin
        ntuple(_->Matrix{ComplexF64}(undef, size(us, 1), size(us, 1)), Val(D))
    end
) where {D}
    Es, us = solve(ptbm, k; bloch_phase=Val(false))
    Eₙ = Es[n]
    if count(≈(Eₙ), Es) > 1
        error(lazy"band $n is degenerate at $k: non-Abelian Berry curvature is not presently handled")
    end
    ∇ptbm = gradient_wrt_momentum(ptbm)

    # NB: Below we assume and require H to be Hermitian
    if D == 2
        return _berrycurvature_2d(k, ∇ptbm, Es, us, n, ∇Hs)
    elseif D == 3
        return _berrycurvature_3d(k, ∇ptbm, Es, us, n, ∇Hs)
    else
        error(lazy"Berry curvature is not implemented for $D-dimensional systems")
    end
end

function _berrycurvature_2d!(k, ∇ptbm, Es, us, n, ∇Hs::NTuple{2, Matrix{ComplexF64}})
    components = (1, 2) # k₁ & k₂ components of Hamiltonian gradient
    ∇Hs = ∇ptbm(k, components, ∇Hs)

    Eₙ = Es[n]
    uₙ = @view us[:, n]
    Ω³ = 0.0
    for m in eachindex(Es)
        m == n && continue
        uₘ = @view us[:, m]
        Eₘ = Es[m]
        Ω³ -= 2*imag(dot(uₙ, ∇Hs[1], uₘ) * dot(uₘ, ∇Hs[2], uₙ)) / (Eₙ - Eₘ)^2  # Ω³ = Ω₁₂
    end

    return Ω³
end

function _berrycurvature_3d!(k, ∇ptbm, Es, us, n, ∇Hs::NTuple{3, Matrix{ComplexF64}})
    components = (1, 2, 3) # k₁, k₂, & k₃ components of Hamiltonian gradient
    ∇Hs = ∇ptbm(k, components, ∇Hs)

    Eₙ = Es[n]
    uₙ = @view us[:, n]
    Ω¹ = Ω² = Ω³ = 0.0
    for m in eachindex(Es)
        m == n && continue
        uₘ = @view us[:, m]
        Eₘ = Es[m]

        Ω¹ -= 2*imag(dot(uₙ, ∇Hs[2], uₘ) * dot(uₘ, ∇Hs[3], uₙ)) / (Eₙ - Eₘ)^2  # Ω¹ = Ω₂₃
        Ω² -= 2*imag(dot(uₙ, ∇Hs[3], uₘ) * dot(uₘ, ∇Hs[1], uₙ)) / (Eₙ - Eₘ)^2  # Ω² = Ω₃₁
        Ω³ -= 2*imag(dot(uₙ, ∇Hs[1], uₘ) * dot(uₘ, ∇Hs[2], uₙ)) / (Eₙ - Eₘ)^2  # Ω³ = Ω₁₂
    end

    return SVector{3, Float64}(Ω¹, Ω², Ω³)
end