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

# Implementation

The Berry curvature is evaluated using the Kubo-like or "sum-over-states" formula [^1]:

```math
\\Omega_{ij}^n
=
\\sum_{n' \\neq n}
\\frac{\\langle n |\\partial H/\\partial k_i|n'\\rangle \\langle n'|\\partial H/\\partial k_j|n\\rangle
       - (j \\leftrightarrow i)}
      {(E_n - E_{n'})^2}
=
-2\\mathrm{Im} \\sum_{n' \\neq n}
\\frac{\\langle n |\\partial H/\\partial k_i|n'\\rangle \\langle n'|\\partial H/\\partial k_j|n\\rangle}
      {(E_n - E_{n'})^2},
```

with ``\\Omega_{ij}^n`` denoting  the Berry curvature tensor component of the `n`th
band and ``E_n`` denoting the energy of the `n`th band. The Berry curvature pseudovector and
(antisymmetric) tensor components are related by ``\\Omega_1 = \\Omega_{23}``,
``\\Omega_2 = \\Omega_{31}``, and ``\\Omega_3 = \\Omega_{12}``.

!!! note
    The sum-over-states formula is of course equivalent to the more commmonly found
    definition involving derivatives of the Berry connection. In particular, the above
    expression corresponds to a definition of the Berry curvature tensor as:

    ```math
    \\Omega_{ij}^n = \\partial_{k_i} A_j^n - \\partial_{k_j} A_i^n,
    ```

    with Berry connection ``A_i^n = \\mathrm{i} \\langle u_{n\\mathbf{k}} | 
    \\nabla_{\\mathbf{k}} | u_{n\\mathbf{k}} \\rangle``. In terms of the Berry curvature
    pseudovector, this is equivalent to ``\\boldsymbol{\\Omega}^n = \\nabla_{\\mathbf{k}}
    \\times \\mathrm{i} \\langle u_{n\\mathbf{k}} | \\nabla_{\\mathbf{k}} |
    u_{n\\mathbf{k}} \\rangle``.

    Other sign conventions exist in the literature (e.g., opposite sign), but the above
    is prevailing in the physics and condensed-matter literature.

[^1]: Xiao, Chang, & Niu, *Berry phase effects on electronic properties*, 
    [Rev. Mod. Phys. **82**, 1959 (2010)](https://doi.org/10.1103/RevModPhys.82.1959), 
    Eq. (1.13).

# Example

Consider the following model in plane group **p**2, without time-reversal symmetry:
```jldoctest berry-curvature
julia> brs = calc_bandreps(2, Val(2); timereversal=false);

julia> cbr = @composite brs[1] + brs[3]
8-irrep CompositeBandRep{2}:
 (1d|A) + (1c|A) (2 bands)

julia> tbm = tb_hamiltonian(cbr, [[0,0], [1,0]]);
```
`tbm` is an 8-term model, whose 6th and 8th terms break time-reversal symmetry, but which
retains inversion symmetry. We can build a parameterized instance of this model and compute
its Berry curvature at some **k**-point to verify that its evenness under inversion (`k` →
`-k`):
```jldoctest berry-curvature; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> ptbm = tbm([0, .1, -0, 0, .2, 0, 0, .3]);

julia> k = [0.1, 0.2];

julia> C₁_k₊ = berrycurvature(ptbm, k, 1)  # Berry curvature of band 1 at `k`
-0.2529580943483421

julia> C₁_k₋ = berrycurvature(ptbm, -k, 1) # Berry curvature of band 1 at `-k`
-0.2529580943483421

julia> C₁_k₊ ≈ C₁_k₋
true

julia> C₁_k₊ ≈ -berrycurvature(ptbm, k, 2) # opposite Berry curvatures in bands 1 and 2
true
```

The model also happens to have a nonzero Chern number, as can be verified by integrating
the Berry curvature over the Brillouin zone using [`chern`](@ref), or manually by summation:
```jldoctest berry-curvature; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> chern(ptbm, 1, 51)
-1.000000071894951
```

# Extended help

As elsewhere in SymmetricTightBinding.jl, this function operates in *reduced* coordinates,
with the momentum `k = (k₁, k₂, k₃)` indicating coefficients relative to the primitive
reciprocal lattice vectors ``\\mathbf{b}_i``, i.e., ``\\mathbf{k} = k_1\\mathbf{b}_1 +
k_2\\mathbf{b}_2 + \\ldots``.
Importantly, this also implies that returned Berry curvature components are *not* the
Cartesian components, but components relative to the primitive reciprocal basis.

**In 3D**, this implies that the returned vector `Ω = [Ω₁, Ω₂, Ω₃]` relates to the Berry
curvature pseudovector ``\\boldsymbol{\\Omega}`` via

```math
\\boldsymbol{\\Omega}
=
V_{\\text{BZ}}^{-1} (\\Omega_1 \\mathbf{b}_1 + \\Omega_2 \\mathbf{b}_2 + \\Omega_3 \\mathbf{b}_3)
=
V_{\\text{BZ}}^{-1} (\\Omega_{23} \\mathbf{b}_1 + \\Omega_{31} \\mathbf{b}_2 + \\Omega_{12} \\mathbf{b}_3),
```

with ``V_{\\text{BZ}} = \\mathbf{b}_1 \\cdot (\\mathbf{b}_2 \\times \\mathbf{b}_3)``
denoting the (signed) volume of the 3D Brillouin zone.
Here, the tensor components ``\\Omega_{12}`` etc. correspond to momentum derivatives in the
*reduced* coordinates, i.e., the derivatives ``\\partial/\\partial k_{1,2}`` are taken with
respect to  reduced (dimensionless) coordinates ``k_{1,2}`` rather than Cartesian ones.

As a result, the returned Berry curvature pseudovector components ``[Ω₁, Ω₂, Ω₃]`` are
related to the Cartesian pseudovector components ``[Ω_x, Ω_y, Ω_z]`` by the primitive
reciprocal lattice matrix ``\\mathbf{B} = [\\mathbf{b}_1 \\mathbf{b}_2 \\mathbf{b}_3]``:

```math
[Ω_x, Ω_y, Ω_z]^{\\text{T}} = V_{\\text{BZ}}^{-1} \\mathbf{B} [Ω₁, Ω₂, Ω₃]^{\\text{T}}
```

Equivalently, the antisymmetric Berry curvature tensors ``\\Omega_{ij}`` (``i,j \\in
\\{1, 2, 3\\}``, i.e., reduced coordinate basis) and ``\\Omega_{\\mu\\nu} (``\\mu, \\nu = 
\\{x, y, z\\}``, i.e. Cartesian basis) relates via

```math
\\Omega_{ij} = \\sum_{\\mu\\nu} B^T_{i\\mu} \\Omega_{\\mu\\nu} B_{\\nu j}\\
\\Omega_{\\mu\\nu} = \\sum_{ij} B^{-T}_{\\mu i} \\Omega_{ij} B^{-1}_{j\\nu}.
```

**In 2D**, the returned Berry curvature is a (pseudo)scalar ``Ω₃ = Ω₁₂``. Again, this value
represents the component of the Berry curvature in the reduced coordinate system. It
is related to the Cartesian pseudoscalar ``\\Omega_{xy}`` by the (signed) area of the
Brillouin zone ``A_{\\text{BZ}} = (\\mathbf{b}_1 \\times \\mathbf{b}_2) \\cdot 
\\hat{\\mathbf{z}}``:

```math
\\Omega_z = \\Omega_{xy} = \\Omega_{12} / A_{\\text{BZ}}.
```

### Implications for surface integrals and flux

When using the returned 3D Berry curvature pseudovector `[Ω₁, Ω₂, Ω₃]` to calculate
physical observables like the Chern number of a slice or the flux through a surface, it's
crucial to integrate the correct physical quantity. The assumed goal is always to compute
a flux integral:

```math
\\Phi = \\int_S \\boldsymbol{\\Omega} \\cdot \\,\\mathrm{d}\\mathbf{S}.
```

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
V_{\\text{BZ}}^{-1} 
(\\Omega_1 \\mathbf{b}_1 + \\Omega_2 \\mathbf{b}_2 + \\Omega_3 \\mathbf{b}_3)
\\cdot
(\\mathbf{b}_1 \\times \\mathbf{b}_2) \\, \\mathrm{d}k_1 \\mathrm{d}k_2.
```

By the properties of the triple product, the terms proportional to ``\\Omega_{1,2}`` vanish,
leaving:

```math
\\mathrm{d}\\Phi
=
V_{\\text{BZ}}^{-1} \\Omega_3 (\\mathbf{b}_3 \\cdot (\\mathbf{b}_1 \\times \\mathbf{b}_2))
\\, \\mathrm{d}k_1 \\mathrm{d}k_2.
```

The term ``\\mathbf{b}_3 \\cdot (\\mathbf{b}_1 \\times \\mathbf{b}_2)`` is the volume of the
3D Brillouin zone, ``V_{\\text{BZ}}``, cancelling the ``V_{\\text{BZ}}^{-1}`` prefactor.
Thus, the flux integral simplifies to:

```math
\\mathrm{d}\\Phi = \\Omega_3 \\, \\mathrm{d}k_1 \\mathrm{d}k_2.
```

In other words: while the distinction between e.g., ``\\Omega_{xy}`` and ``\\Omega_{12}``
is important for interpreting the curvature at an individual **k**-point, it can be ignored
when calculating the Chern number as an integral over the reduced coordinates (e.g., 
``k_{1,2} \\in (-1/2, 1/2]`` above), since the transformation factors in the curvature
and the area element ``\\mathrm{d}\\mathbf{S}`` cancel out.
"""
function berrycurvature(
    ptbm::ParameterizedTightBindingModel{D},
    k::ReciprocalPointLike{D},
    n::Integer,
    ∇Hs::NTuple{D, Matrix{ComplexF64}} = begin
        ntuple(_ -> Matrix{ComplexF64}(undef, ptbm.tbm.N, ptbm.tbm.N), Val(D))
    end
) where {D}
    1 ≤ n ≤ ptbm.tbm.N || error("band index `n` is beyond the model's size")
    Es, us = solve(ptbm, k; bloch_phase=Val(false))
    Eₙ = Es[n]
    if count(≈(Eₙ), Es) > 1
        error(lazy"band $n is degenerate at $k: non-Abelian Berry curvature is not presently handled")
    end
    ∇ptbm = gradient_wrt_momentum(ptbm)

    # NB: Below we assume and require H to be Hermitian
    if D == 2
        return _berrycurvature_2d!(k, ∇ptbm, Es, us, n, ∇Hs)
    elseif D == 3
        return _berrycurvature_3d!(k, ∇ptbm, Es, us, n, ∇Hs)
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


"""
    chern(
        ptbm::ptbm::ParameterizedTightBindingModel{2},
        n::Integer,
        Nk::Integer,
    ) --> Float64

Evaluate the Chern number of the `n`th band of the 2D model `ptbm`, by grid-based
integration of the Berry curvature (using [`berrycurvature`](@ref)) over the
reduced-coordinate Brillouin zone ``(-1/2, 1/2]×(-1/2, 1/2]`` using `Nk`×`Nk` sample points.

Note that, unlike the Fukui et al. method ([`chern_fukui`](@ref)), this approach does not
guarantee integer results for finite `Nk`, but will converge smoothly to the correct integer
Chern number with increasing `Nk`.

## Definition

The returned Chern number is defined with the sign-convention:

```math
C = \\frac{1}{2\\pi} \\int_{\\text{BZ}} \\nabla_{\\mathbf{k}} \\times \\mathrm{i}
    \\langle u_{n\\mathbf{k}} | \\nabla_{\\mathbf{k}} | u_{n\\mathbf{k}} \\rangle
    \\, \\mathrm{d}^2 \\mathbf{k}.
```
"""
function chern(ptbm::ParameterizedTightBindingModel{2}, n::Integer, Nk::Integer)
    ks = range(-0.5, 0.5, Nk+1)[2:end]
    ∇Hs = ntuple(_ -> Matrix{ComplexF64}(undef, ptbm.tbm.N, ptbm.tbm.N), Val(2))
    Ωs = (berrycurvature(ptbm, (k1,k2), n, ∇Hs) for k1 in ks, k2 in ks) # berry curvatures
    Φ = sum(Ωs; init = zero(Float64)) / Nk^2 # flux
    return Φ / (2π) # Chern number
end

"""
    chern_fukui(
        ptbm::ptbm::ParameterizedTightBindingModel{2},
        n::Integer,
        Nk::Integer,
    ) --> Int

Evaluate the Chern number of the `n`th band of the 2D model `ptbm`, using the link-variable
integration approach of Fukui et al. [^1] and `Nk` linear sample points per **k**-space
dimension.

The returned value is guaranteed to be an integer (but will only converge to the correct
Chern number for sufficiently high `Nk`).

[^1]: Fukui, Hatsugai, & Suzuki, *Chern Numbers in Discretized Brillouin Zone: Efficient
     Method of Computing (Spin) Hall Conductances*,
     [J. Phys. Soc. Jpn. **74**, 1674 (2005)](https://doi.org/10.1143/JPSJ.74.1674).
     
See also [`chern`](@ref) for an alternative approach based on direct integration of the
Berry curvature, evaluated from a sum-over-states formula. The Fukui approach is usually
more performant for Chern number evaluation.

## Definition

The returned Chern number is defined with the sign-convention:

```math
C = \\frac{1}{2\\pi} \\int_{\\text{BZ}} \\nabla_{\\mathbf{k}} \\times \\mathrm{i}
    \\langle u_{n\\mathbf{k}} | \\nabla_{\\mathbf{k}} | u_{n\\mathbf{k}} \\rangle
    \\, \\mathrm{d}^2 \\mathbf{k}.
```
"""
function chern_fukui(ptbm::ParameterizedTightBindingModel{2}, n::Integer, Nk::Integer)
    # single-band (Abelian) Fukui approach
    1 ≤ n ≤ ptbm.tbm.N || error("band index `n` is beyond the model's size")
    Nk == 1 && error("`Nk` must be at least 2 to define a meaningful Fukui discretization")
    ks = range(-0.5, 0.5, Nk+1)
    # NB: ↑ We include the endpoint here, since we do _not_ have H(k) = H(k+G) in Convention
    #     1: trying to naively reuse the eigenvectors at kⱼ = -0.5 for kⱼ = 0.5 would be
    #     wrong in this convention (since, in this convention they differ by a gauge phase
    #     of exp(iG·qᵢ)): while we could exploit this to avoid a "boundary solve", it would
    #     specialize the implementation to Convention 1, failing in Convention 2. By just
    #     not wrapping the k-values explicitly across the BZ, the implementation works in
    #     both conventions.

    # Algorithm & implementation notes:
    # Below is the standard Fukui algorithm, running over plaquettes with corners at
    # (kᵢ, kⱼ), (kᵢ₊₁, kⱼ), (kᵢ, kⱼ₊₁), (kᵢ₊₁, kⱼ₊₁); the code is a bit more involved than
    # a naive implementation, only to avoid redundant eigen-solves; we do this by storing
    # two rows of eigenvectors at a time (for fixed kⱼ and kⱼ₊₁), and reusing the `j`th row
    # when moving to the next `j+1`th row; this ensures we only do (Nk+1)² solves,
    # rather than 4(Nk+1)² solves (the reason it is not Nk², is that we do not exploit the
    # wrap-around in either the i or j directions; exploiting it in j would also require
    # an additional row-cache, which is just too much hassle for asymptotically little gain)

    eigen_storeⱼ = Matrix{ComplexF64}(undef, ptbm.tbm.N, length(ks)) # eigvecs of the `j`th row
    k₁ = ks[1]
    for i in eachindex(ks) # initialize first row at `j = 1`
        kᵢ = ks[i]
        vs = solve(ptbm, (kᵢ, k₁))[2]
        eigen_storeⱼ[:, i] .= @view vs[:, n]
    end

    eigen_storeⱼ₊₁ = similar(eigen_storeⱼ) # eigvecs of the `j+1`th row
    C = zero(Float64)
    for j in Base.OneTo(length(ks)-1)
        jp1 = j+1
        kⱼ₊₁ = ks[jp1]
        
        for i in eachindex(ks) # compute next row of eigvecs at `j+1`
            kᵢ = ks[i]
            vs = solve(ptbm, (kᵢ, kⱼ₊₁))[2]
            eigen_storeⱼ₊₁[:,i] .= @view vs[:, n]
        end

        # now compute flux in current plaquette, bounded by (i, j, i+1, j+1)
        for i in Base.OneTo(length(ks)-1)
            ip1 = i+1
            
            v₀₀ = @view eigen_storeⱼ[:,i]     # (kᵢ, kⱼ)
            v₁₀ = @view eigen_storeⱼ[:,ip1]   # (kᵢ₊₁, kⱼ)
            v₀₁ = @view eigen_storeⱼ₊₁[:,i]   # (kᵢ, kⱼ₊₁)
            v₁₁ = @view eigen_storeⱼ₊₁[:,ip1] # (kᵢ₊₁, kⱼ₊₁)
            
            # NB: Note the sign convention below: Fukui et al. (and mathematicians, broadly)
            #     define the Berry connection to be Aᴵ = ⟨uₖ|∇ₖ|uₖ⟩ (an imaginary quantity),
            #     but the usual physicist definition is Aᴵᴵ = i⟨uₖ|∇ₖ|uₖ⟩ (a real quantity).
            #     Similarly, Fukui defines the Chern number as Cᴵ = (2πi)⁻¹ ∫ ∇ₖ×Aᴵ d²k, 
            #     while physicists generally use Cᴵᴵ = (2π)⁻¹ ∫ ∇ₖ×Aᴵᴵ d²k. Taken together,
            #     there's a difference Cᴵᴵ = -Cᴵ. We prefer the physicist convention and
            #     take Ω = ∇ₖ×Aᴵᴵ as the Berry curvature and C = (2π)⁻¹ ∫ Ω d²k = Cᴵᴵ as the
            #     Chern number. This just means that we need to negate Fukui et al.'s
            #     expressions.
            U_□ = dot(v₀₀, v₁₀) * dot(v₁₀, v₁₁) * dot(v₁₁, v₀₁) * dot(v₀₁, v₀₀)
            C -= atan(imag(U_□), real(U_□))
            # NB: the above expression is equivalent to `C -= imag(log(z))`, since 
            #     log(z) = log|z| + i arg(z), so we can use `atan(imag(U_□), real(U_□))`
            #     (which gives arg(U_□)) instead to be slightly more efficient
        end
        # swap `j`th and `j+1`th rows for next iteration (`j+1`th row becomes new `j`th row)
        eigen_storeⱼ, eigen_storeⱼ₊₁ = eigen_storeⱼ₊₁, eigen_storeⱼ
    end

    C /= 2π

    # C is a floating point number now, but is guaranteed to be numerically close to an
    # integer (only deviating due to numerical imprecision (in e.g., `log` and eigen-solves,
    # but not from k-space discretization): before returning, we round to nearest integer
    Cint = round(Int, C)
    abs(C - Cint) > 1e-4 && error(lazy"computed Chern number $C deviates more substantially from an integer than expected by the implementation")
    return Cint
end


function chern_fukui(
    ptbm::ParameterizedTightBindingModel{2},
    nsv::AbstractVector{T},
    Nk::Integer
) where T<:Union{Integer, AbstractRange{<:Integer}}
    # exactly the same idea as single-band (Abelian) Fukui approach, but now for non-Abelian
    # case and possibly with multiple bands

    Nk == 1 && error("`Nk` must be at least 2 to define a meaningful Fukui discretization")
    ks = range(-0.5, 0.5, Nk+1)

    nsv_flat = Iterators.flatten(nsv)
    allunique(nsv_flat) || error("requested band multiplets `nsv` have overlapping band indices")
    Nb = sum(length, nsv_flat; init = 0)

    eigen_storeⱼ = Array{ComplexF64}(undef, ptbm.tbm.N, Nb, length(ks)) # eigvecs of the `j`th row
    k₁ = ks[1]
    for i in eachindex(ks) # initialize first row at `j = 1`
        kᵢ = ks[i]
        vs = solve(ptbm, (kᵢ, k₁))[2]
        for (b, n) in enumerate(nsv_flat)
            eigen_storeⱼ[:, b, i] .= @view vs[:, n]
        end
    end

    eigen_storeⱼ₊₁ = similar(eigen_storeⱼ) # eigvecs of the `j+1`th row
    Cs = zeros(Float64, length(nsv))
    for j in Base.OneTo(length(ks)-1)
        jp1 = j+1
        kⱼ₊₁ = ks[jp1]
        
        for i in eachindex(ks) # compute next row of eigvecs at `j+1`
            kᵢ = ks[i]
            vs = solve(ptbm, (kᵢ, kⱼ₊₁))[2]
            for (b, n) in enumerate(nsv_flat)
                eigen_storeⱼ₊₁[:, b, i] .= @view vs[:, n]
            end
        end

        for i in Base.OneTo(length(ks)-1) # computes fluxes in plaquette (i, j, i+1, j+1)
            ip1 = i+1
            
            b = 1
            for (r, ns) in enumerate(nsv)
                N_r = length(ns)
                if N_r == 1 # scalar case (Abelian)
                    v₀₀ = @view eigen_storeⱼ[:, b, i]     # (kᵢ, kⱼ)
                    v₁₀ = @view eigen_storeⱼ[:, b, ip1]   # (kᵢ₊₁, kⱼ)
                    v₀₁ = @view eigen_storeⱼ₊₁[:, b, i]   # (kᵢ, kⱼ₊₁)
                    v₁₁ = @view eigen_storeⱼ₊₁[:, b, ip1] # (kᵢ₊₁, kⱼ₊₁)
                    U_□ = dot(v₀₀, v₁₀) * dot(v₁₀, v₁₁) * dot(v₁₁, v₀₁) * dot(v₀₁, v₀₀)
                    Cs[r] -= atan(imag(U_□), real(U_□))

                else # matrix case (non-Abelian)
                    bs = b:b+N_r-1
                    vs₀₀ = @view eigen_storeⱼ[:, bs, i]     # (kᵢ, kⱼ)
                    vs₁₀ = @view eigen_storeⱼ[:, bs, ip1]   # (kᵢ₊₁, kⱼ)
                    vs₀₁ = @view eigen_storeⱼ₊₁[:, bs, i]   # (kᵢ, kⱼ₊₁)
                    vs₁₁ = @view eigen_storeⱼ₊₁[:, bs, ip1] # (kᵢ₊₁, kⱼ₊₁)
                    U_□ = (vs₀₀' * vs₁₀) * (vs₁₀' * vs₁₁) * (vs₁₁' * vs₀₁) * (vs₀₁' * vs₀₀)
                    Cs[r] -= imag(tr(log(U_□))) # tr(log(…) for gauge-invariant contribution
                end
                b += N_r
            end
        end
        eigen_storeⱼ, eigen_storeⱼ₊₁ = eigen_storeⱼ₊₁, eigen_storeⱼ
    end

    Cs /= 2π

    return round.(Int, Cs)
end