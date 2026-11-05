"""
    quantummetric(
        ptbm::ParameterizedTightBindingModel{D},
        k::ReciprocalPointLike{D},
        ns::Union{Integer, AbstractUnitRange{<:Integer}},
        [∇Hs::NTuple{D, Matrix{ComplexF64}}]
    )  --> g :: Hermitian{Float64, <:SMatrix{D, D, Float64}}

Compute the quantum metric `g` of the `ns`th band - or, if `ns` is a contiguous range of band
indices, of the band multiplet `ns` (see "Band multiplets" below) - of the
coefficient-parameterized tight-binding model `ptbm` at the **k**-point `k`. The model must
be Hermitian.
The returned `g` is a real, symmetric, positive semi-definite `D`×`D` matrix.

The optional `∇Hs` argument is a tuple of work-matrices used to store the momentum gradient
of the Hamiltonian.

# Implementation

The quantum metric is evaluated using the Kubo-like or "sum-over-states" formula [^1]:

```math
g_{ij}^n
=
\\mathrm{Re} \\sum_{n' \\neq n}
\\frac{\\langle n |\\partial H/\\partial k_i|n'\\rangle \\langle n'|\\partial H/\\partial k_j|n\\rangle}
      {(E_n - E_{n'})^2},
```

with ``g_{ij}^n`` denoting the quantum metric of the `n`th band and ``E_n`` its energy.
This is the same expression used to compute the Berry curvature (see
[`berrycurvature`](@ref)), but taking the real rather than the imaginary part: both are
parts of a single object, the quantum geometric tensor
``Q_{ij} = g_{ij} - \\tfrac{\\mathrm{i}}{2}\\Omega_{ij}`` (see
[`quantumgeometrictensor`](@ref)).

[^1]: Provost & Vallee, *Riemannian structure on manifolds of quantum states*,
    [Commun. Math. Phys. **76**, 289 (1980)](https://doi.org/10.1007/BF02193559).

## Band multiplets

If `ns` is a contiguous range of band indices ``\\mathcal{N} = \\{n_1, \\ldots, n_M\\}``,
the sum over intermediate states is restricted to bands *outside* ``\\mathcal{N}``:

```math
g_{ij}^{\\mathcal{N}}
=
\\mathrm{Re} \\sum_{n \\in \\mathcal{N}} \\sum_{m \\notin \\mathcal{N}}
\\frac{\\langle n |\\partial H/\\partial k_i|m\\rangle \\langle m|\\partial H/\\partial k_j|n\\rangle}
      {(E_n - E_m)^2}
=
\\tfrac{1}{2}\\mathrm{Tr}[(\\partial_{k_i} P)(\\partial_{k_j} P)],
```

with ``P = \\sum_{n \\in \\mathcal{N}} |u_{n\\mathbf{k}}\\rangle\\langle
u_{n\\mathbf{k}}|`` the projector onto the multiplet. The right-hand form, built from `P`
alone, shows that this is the part of the geometry unchanged by any mixing of the bands in
``\\mathcal{N}`` among themselves. The restriction to ``m \\notin \\mathcal{N}`` is not an
ad-hoc omission: it is what the projector ``1 - P`` produces, and it is also why the
intra-multiplet - i.e., the gauge-dependent and degeneracy-divergent - terms never arise.

Unlike the Berry curvature, the metric is *not* additive over the multiplet,
``g^{\\mathcal{N}} \\neq \\sum_{n \\in \\mathcal{N}} g^n``: the omitted intra-multiplet terms
cancel pairwise under the imaginary part taken by [`berrycurvature`](@ref), but add under the
real part taken here.

Since ``g^{\\mathcal{N}}`` depends only on ``P``, degeneracies *within* `ns` are
unproblematic, and the multiplet form is the only well-defined choice where bands touch; a
degeneracy between a band in `ns` and a band outside `ns` is an error.

# Example

```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> brs = bandreps(2, Val(2); timereversal=false);

julia> cbr = @composite brs[1] + brs[3]; # (1d|A) + (1c|A)

julia> tbm = tb_hamiltonian(cbr, [[0,0], [1,0]]);

julia> ptbm = tbm([0, .1, 0, 0, .2, 0, 0, .3]);

julia> quantummetric(ptbm, [0.1, 0.2], 1)
2×2 Hermitian{Float64, SMatrix{2, 2, Float64, 4}} with indices SOneTo(2)×SOneTo(2):
  0.0713504   -0.00803914
 -0.00803914   0.225108
```

# Extended help

As for [`berrycurvature`](@ref), the returned components refer to *reduced* coordinates,
i.e., ``g_{ij}`` involves derivatives with respect to the reduced coordinates ``k_{1,2,3}``
rather than Cartesian ones.
Conversion from a reduced setting, ``g_{ij}``, to a Cartesian one, ``g_{\\mu\\nu}``,
parallels the procedure described in [`berrycurvature`](@ref):
``g_{\\mu\\nu} = \\sum_{ij} B^{-\\mathrm{T}}_{\\mu i} g_{ij} B^{-1}_{j\\nu}`` with
``\\mathbf{B} = [\\mathbf{b}_1 \\ldots \\mathbf{b}_D]``.

Unlike the Chern number, however, this conversion cannot generally be omitted: there, the
transformation factors of ``\\Omega`` and of the area element cancel, but no such cancellation
occurs for the metric. In particular, contractions such as ``\\mathrm{tr}\\, g`` are
coordinate-dependent: ``\\sum_i g_{ii}`` in reduced coordinates is not rotation-invariant,
whereas the Cartesian trace is, equaling ``\\sum_{ij} (\\mathbf{G}^{-1})_{ij} g_{ij}`` with
``G_{ij} = \\mathbf{b}_i \\cdot \\mathbf{b}_j`` the metric tensor of the reciprocal lattice.

!!! note "Dependence on the orbital embedding"
    The quantum metric depends on where the orbitals sit in the unit cell - unlike the Chern
    number, which does not. The returned values are those of the Convention 1 Fourier
    transform used throughout SymmetricTightBinding, i.e., they include the orbital-position
    phases set by the Wyckoff positions of the underlying band representation; Convention 2
    values will in general differ (see [Fourier transforms conventions](@ref)).
"""
function quantummetric(
    ptbm::ParameterizedTightBindingModel{D},
    k::ReciprocalPointLike{D},
    ns::Union{Integer, AbstractUnitRange{<:Integer}},
    ∇Hs::NTuple{D, Matrix{ComplexF64}} = begin
        ntuple(_ -> Matrix{ComplexF64}(undef, ptbm.tbm.N, ptbm.tbm.N), Val(D))
    end
) where {D}
    Q = quantumgeometrictensor(ptbm, k, ns, ∇Hs)
    return real(Q) # quantum metric from the quantum geometric tensor: gᵢⱼ = Re Qᵢⱼ
end

"""
    quantumgeometrictensor(
        ptbm::ParameterizedTightBindingModel{D},
        k::ReciprocalPointLike{D},
        ns::Union{Integer, AbstractUnitRange{<:Integer}},
        [∇Hs::NTuple{D, Matrix{ComplexF64}}]
    )  --> Q :: Hermitian{ComplexF64, <:SMatrix{D, D, ComplexF64}}

Compute the quantum geometric tensor `Q` of band `ns` - or, if `ns` is a contiguous range of
band indices, of the band multiplet `ns` - of the coefficient-parameterized tight-binding
model `ptbm` at the **k**-point `k`. The returned `Q` is a complex, Hermitian, positive
semi-definite `D`×`D` matrix, evaluated via the sum-over-states formula:

```math
Q_{ij}^{\\mathcal{N}}
=
\\sum_{n \\in \\mathcal{N}} \\sum_{m \\notin \\mathcal{N}}
\\frac{\\langle n |\\partial H/\\partial k_i|m\\rangle \\langle m|\\partial H/\\partial k_j|n\\rangle}
      {(E_n - E_m)^2}
=
g_{ij}^{\\mathcal{N}} - \\tfrac{\\mathrm{i}}{2}\\Omega_{ij}^{\\mathcal{N}},
```

with ``\\mathcal{N}`` the bands in `ns`. Equivalently, ``Q_{ij}^{\\mathcal{N}}
= \\mathrm{Tr}[(\\partial_{k_i} P)(1 - P)(\\partial_{k_j} P)]`` [^1],
with ``P = \\sum_{n \\in \\mathcal{N}} |u_{n\\mathbf{k}}\\rangle\\langle u_{n\\mathbf{k}}|``
the projector onto the multiplet. That is, `Q` combines the quantum metric ``g_{ij}`` and
the Berry curvature tensor ``\\Omega_{ij}``, as returned by [`quantummetric`](@ref) and
[`berrycurvature`](@ref); see either for details on input arguments, band multiplets, and
coordinate conventions.

[^1]: See e.g., Ma, Chen, Fan, & Liu, *Abelian and non-Abelian quantum geometric tensor*,
    [Phys. Rev. B **81**, 245129 (2010)](https://doi.org/10.1103/PhysRevB.81.245129), Eq.
    (9), for the multiplet formulation, and Graf & Piéchon, *Berry curvature and quantum
    metric in N-band systems: An eigenprojector approach*,
    [Phys. Rev. B **104**, 085114 (2021)](https://doi.org/10.1103/PhysRevB.104.085114),
    Eq. (3) (there, for a single band).
"""
function quantumgeometrictensor(
    ptbm::ParameterizedTightBindingModel{D},
    k::ReciprocalPointLike{D},
    ns′::Union{Integer, AbstractUnitRange{<:Integer}},
    ∇Hs::NTuple{D, Matrix{ComplexF64}} = begin
        ntuple(_ -> Matrix{ComplexF64}(undef, ptbm.tbm.N, ptbm.tbm.N), Val(D))
    end
) where {D}
    S = hermiticity(ptbm)
    S === HERMITIAN || error(lazy"the quantum geometry requires a Hermitian model; got $S")
    ns = ns′ isa Integer ? (ns′:ns′) : ns′
    isempty(ns) && error("band specification `ns` must be non-empty")
    1 ≤ first(ns) && last(ns) ≤ ptbm.tbm.N || error("band index is beyond the model's size")
    Es, us = solve(ptbm, k; bloch_phase=Val(false))
    _check_multiplet_isolated(Es, ns, k)

    ∇ptbm = gradient_wrt_momentum(ptbm)
    components = ntuple(identity, Val(D)) # k₁, …, k_D components of Hamiltonian gradient
    ∇Hs = ∇ptbm(k, components, ∇Hs)

    # sum-over-states formula, wrt. reduced coordinates:
    #       Qᵢⱼ = ∑_{n ∈ ns} ∑_{m ∉ ns} ⟨n|∂H/∂kᵢ|m⟩⟨m|∂H/∂kⱼ|n⟩ / (Eₙ - Eₘ)²
    Q = zero(MMatrix{D, D, ComplexF64})
    for n in ns
        Eₙ = Es[n]
        uₙ = @view us[:, n]
        for m in eachindex(Es)
            m ∈ ns && continue # intermediate states only from outside the multiplet
            uₘ = @view us[:, m]
            Eₘ = Es[m]
            # vᵢ ≡ ⟨n|∂H/∂kᵢ|m⟩; since H is Hermitian, the "partner" matrix element of the
            # sum-over-states numerator is ⟨m|∂H/∂kⱼ|n⟩ = conj(vⱼ), so that the numerator
            # ⟨n|∂H/∂kᵢ|m⟩⟨m|∂H/∂kⱼ|n⟩ = vᵢ conj(vⱼ)
            vs = ntuple(i -> dot(uₙ, ∇Hs[i], uₘ), Val(D))
            ΔE⁻² = 1 / (Eₙ - Eₘ)^2
            for j in 1:D
                for i in 1:j # upper triangle only; lower triangle implied by `Hermitian`
                    Q[i, j] += vs[i] * conj(vs[j]) * ΔE⁻²
                end
            end
        end
    end

    return Hermitian(SMatrix(Q), :U) # Qⱼᵢ = Qᵢⱼ*
end

# error if a band in `ns` is degenerate with a band outside `ns`: the sum-over-states
# denominators then diverge, and the geometry of the multiplet `ns` is undefined. The
# degeneracy tolerance, `max(degen_rtol·bandwidth, degen_atol)`, follows `gradients.jl`; its
# absolute part ensures that degeneracies near `E = 0` are also caught.
function _check_multiplet_isolated(Es, ns, k; degen_rtol = √eps(), degen_atol = 1e-12)
    degen_tol = max(degen_rtol * (maximum(Es) - minimum(Es)), degen_atol)
    for n in ns, m in eachindex(Es)
        (m ∈ ns || abs(Es[m] - Es[n]) ≥ degen_tol) && continue
        error("band $n is degenerate with band $m at $k, but only band(s) $ns were \
               requested: the quantum geometry of a band multiplet is defined only for the \
               multiplet as a whole (extend the band specification to include band $m)")
    end
end

# TODO: natural follow-ups, deliberately left out so far:
#       - Brillouin zone integrals of `tr g`: the Souza-Wilkens-Martin localization marker
#         and the gauge-invariant part Ωᴵ of the Marzari-Vanderbilt spread, in analogy with
#         `chern` (requires a Cartesian contraction)
#       - an overlap-based ("Fukui-like") evaluation from |⟨uₖ|uₖ₊δₖ⟩|², in analogy with
#         `chern_fukui`; more robust than sum-over-states near small gaps
