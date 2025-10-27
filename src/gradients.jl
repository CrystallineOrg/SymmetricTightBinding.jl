struct TightBindingModelHoppingGradient{D}
   tbm :: TightBindingModel{D}
end

"""
    gradient_wrt_hopping(tbm :: TightBindingModel)
    gradient_wrt_hopping(ptbm :: ParameterizedTightBindingModel)

Return a structure that encodes the gradient of a tight-binding model `tbm` or `ptbm` with
respect to the hopping coefficients.

To evaluate the gradient at a particular momentum `k`, use the returned structure as a
functor at `k`. I.e., `gradient(tbm)(k)` returns the gradient of the tight-binding
Hamiltonian with respect to all hoppping coefficients at momentum `k`. This gradient is a
vector of matrices.
"""
gradient_wrt_hopping(tbm::TightBindingModel) = TightBindingModelHoppingGradient(tbm)
gradient_wrt_hopping(ptbm::ParameterizedTightBindingModel) = gradient_wrt_hopping(ptbm.tbm)

function (tbmg::TightBindingModelHoppingGradient{D})(k::ReciprocalPointLike{D}) where D
    map(tbt->tbt(k), tbmg.tbm.terms)
end

# evaluate just the `i`th component of the coefficient gradient
function (tbmg::TightBindingModelHoppingGradient{D})(
    k::ReciprocalPointLike{D},
    i::Int
) where D
    tbmg.tbm[i](k)
end

"""
    energy_gradient_wrt_hopping(
        ptbm::ParameterizedTightBindingModel{D},
        k::ReciprocalPointLike{D}
        (Es, us) = solve(ptbm, k; bloch_phase=Val(false));
        degen_rtol::Float64 = 1e-12,
        degen_atol::Float64 = 1e-12
    ) where D

Return the hopping gradient of the energy of each band in `ptbm` evaluated at momentum `k`.

The gradient is computed using the Feynman-Hellmann theorem. For degenerate bands (assessed
energetically using relative and absolute tolerances `degen_rtol` and `degen_atol`), a
degenerate variant is used, equivalent to degenerate perturbation theory.

The gradient is returned as column vectors, one for each band, with each column containing
the gradient of the corresponding energy with respect to the hopping coefficients of `ptbm`.
"""
function energy_gradient_wrt_hopping(
    ptbm::ParameterizedTightBindingModel{D},
    k::ReciprocalPointLike{D},
    (Es, us) = solve(ptbm, k; bloch_phase=Val(false)) # "unperturbed" energies & eigenstates
    ;
    degen_rtol::Float64 = 1e-12,
    degen_atol::Float64 = 1e-12
) where D
    Nᶜ = length(ptbm.tbm) # number of hopping terms
    Nᵇ = ptbm.tbm.N       # number of bands

    tbmg = gradient_wrt_hopping(ptbm) # ∇ᶜH
    
    # figure out if any bands are degenerate; if so, we need to treat them using degenerate
    # perturbation theory
    degen_tol = max(degen_rtol * (maximum(Es) - minimum(Es)), degen_atol)
    bands = Vector{UnitRange{Int}}(); sizehint!(bands, Nᵇ) # indices into degenerate bands
    n = 1
    while n <= Nᵇ
        if n == Nᵇ # last band
            push!(bands, n:n)
            break
        end
        Eₙ = Es[n]
        n′ = findnext(E->abs(E-Eₙ) < degen_tol, Es, n+1)
        if isnothing(n′) # non-degenerate band
            push!(bands, n:n)
            n += 1
        else             # degenerate bands
            push!(bands, n:n′)
            n = n′ + 1
        end
    end

    # apply Feynman-Hellmann theorem, either in degenerate or non-degenerate variants
    ∇ᶜEs = Matrix{Float64}(undef, Nᶜ, Nᵇ)
    for i in 1:Nᶜ
        ∂ᵢH = tbmg(k, i)
        for ns in bands
            if length(ns) == 1 # non-degenerate band
                n = @inbounds ns[1]
                uₙ = @view us[:, n]
                ∂ᵢEₙ = dot(uₙ, ∂ᵢH, uₙ)
                ∇ᶜEs[i, n] = real(∂ᵢEₙ)
            else               # degenerate bands
                us′ = @view us[:, ns]
                M = us′' * ∂ᵢH * us′ # Mₙₘ = ⟨uₙ|∂ᵢH|uₘ⟩ for n,m ∈ `ns`
                ∂ᵢEs′ = eigvals!(Hermitian(M)) # ∂ᵢEₙ for n in `ns`
                ∇ᶜEs[i, ns] = ∂ᵢEs′
            end
        end
    end

    return eachcol(∇ᶜEs)
end

# ---------------------------------------------------------------------------------------- #
struct TightBindingModelMomentumGradient{D}
   ptbm :: ParameterizedTightBindingModel{D}
end

"""
    gradient_wrt_momentum(ptbm :: ParameterizedTightBindingModel)

Return a structure that encodes the gradient of a parameterized tight-binding model `ptbm`
with respect to its momentum coordinates (in the basis of the primitive reciprocal lattice
vectors)

To evaluate the gradient at a particular momentum `k`, use the returned structure as a
functor at `k`. I.e., `gradient_wrt_momentum(ptbm)(k)[i]` returns the `i`th component of the
momentum derivative of `ptbm` with respect to the momentum at `k`. The return value is a
`D`-dimensional tuple of matrices (see also [`TightBindingModelMomentumGradient`](@ref)).
"""
function gradient_wrt_momentum(ptbm::ParameterizedTightBindingModel{D}) where {D}
    return TightBindingModelMomentumGradient{D}(ptbm)
end

"""
    (∇ptbm::TightBindingModelMomentumGradient{D})(
        k::ReciprocalPointLike{D},
        components::NTuple{C, Int},
        [∇Hs::NTuple{C, Matrix{ComplexF64}}]
    ) --> ∇Hs

Evaluate components of the momentum gradient `∇ptbm` at momentum `k`, writing into the
mutated matrix tuple `∇Hs` (automatically initialized if not provided). 

The momentum components are specified by the `components` argument. These specifications are
assumed relative to the primitive reciprocal lattice vectors ``\\mathbf{b}_i`` (which `k =
(k₁, k₂, k₃)` is also assumed specified relative to, such that the momentum vector is
``\\mathbf{k} = k_1 \\mathbf{b}_1 + k_2 \\mathbf{b}_2 + k_3 \\mathbf{b}_3``).
E.g., for a 3D model with `components = (2, 3)`, the elements of the returned `∇Hs` tuple
have the interpretation that `∇Hs[1]` is the gradient of the Hamiltonian with respect to
`k₂`, i.e., ``\\partial H/\\partial k_2``, and `∇Hs[2]` is the gradient with respect to
`k₃`, ``\\partial H/\\partial k_3``.

## Coordinate transformation

Note that because `k` is specified relative to the primitive reciprocal basis, which is
generally neither orthogonal nor unit-normalized, the computed gradient is not the gradient
with respect to the Cartesian components of ``\\mathbf{k}``, nor with respect to the
coefficients of unit-normalized primitive reciprocal lattice vectors.

One may readily convert to either of these cases by using the chain rule, however. E.g., the
gradient components ``\\partial H/\\partial k̃_i`` for coordinates ``k̃_i`` relative to a
unit-normalized reciprocal basis specification of the momentum vector
``\\mathbf{k} = k̃_1\\hat{\\mathbf{b}}_1 + k̃_2\\hat{\\mathbf{b}}_2 + k̃_3\\hat{\\mathbf{b}}_3``
are related to the returned components by
``\\partial H/\\partial k̃_i = |\\mathbf{b}_i|^{-1} \\partial H/\\partial k_i``.

Similarly, the components relative to Cartesian coordinates ``kᶜ_i`` of the momentum vector
``\\mathbf{k} = kᶜ_1\\hat{\\mathbf{e}}_1 + kᶜ_2\\hat{\\mathbf{e}}_2 + 
kᶜ_3\\hat{\\mathbf{e}}_3`` can be obtained via ``\\partial H/\\partial kᶜ_i = \\sum_j
(\\mathbf{B}^{-\\mathrm{T}})_{ij} \\partial H/\\partial k_j``,
where ``\\mathbf{B} = [\\mathbf{b}_1 \\mathbf{b}_2 \\mathbf{b}_3]`` is a matrix whose
columns are the primitive reciprocal lattice vectors and ``\\mathbf{B}^{-\\mathrm{T}}`` is
its inverse transpose.
"""
function (∇ptbm::TightBindingModelMomentumGradient{D})(
    k::ReciprocalPointLike{D},
    components::NTuple{C, Int},
    ∇Hs::NTuple{C, Matrix{ComplexF64}} = begin
        N = ∇ptbm.ptbm.tbm.N
        ntuple(_ -> zeros(ComplexF64, (N, N)), Val(C))
    end
) where {D, C}
    if length(k) ≠ D
        error("momentum `k` must be a $D-dimensional vector to match the model dimension")
    end
    ptbm = ∇ptbm.ptbm
    tbm = ptbm.tbm
    N = tbm.N
    foreach(∇Hs) do ∇H
        size(∇H) == (N, N) || _throw_scratch_size_mismatch(∇H, N)
        fill!(∇H, zero(ComplexF64)) # make sure storage is reset
    end

    for (tbt, c) in zip(tbm.terms, ptbm.cs) # ↓ modifies `∇Hs` in-place
        evaluate_tight_binding_momentum_gradient_term!(tbt, k, components, c, ∇Hs)
    end

    return ∇Hs
end

"""
    evaluate_tight_binding_momentum_gradient_term!(
        tbt::TightBindingTerm,
        k::ReciprocalPointLike, 
        components::NTuple{N, Int},
        [c::Union{Nothing, <:Number} = nothing],
        [∇Hs::NTuple{N, Matrix{ComplexF64}} = ntuple(zeros(ComplexF64, size(tbt)), N)]
    )

Evaluate components of the momentum gradient of the tight-binding term `tbt` at
momentum `k` over momentum, possibly multiplied by a scalar coefficient `c` (unity if
omitted).
The components are specified by the `components` tuple: a value of `i` for some `components`
element indicates the momentum gradient along `k[i]` (`components` values must be in the
range `1:D`).

The `components[idx]` term is _added_ into the scratch space matrix `∇Hs[idx]`; if `∇Hs`
is not provided, it is initialized as a suitable tuple of zero matrices of appropriate size.
The function returns the modified `∇Hs` matrix tuple.

The function is analogous to `evaluate_tight_binding_term!`, but computes momentum
gradient components rather than the Hamiltonian matrix itself.
"""
function evaluate_tight_binding_momentum_gradient_term!(
    tbt::TightBindingTerm{D},
    k::ReciprocalPointLike{D},
    components::NTuple{C, Int},
    c::Union{Nothing, <:Number} = nothing,
    ∇Hs::NTuple{C, Matrix{ComplexF64}} = ntuple(_ -> zeros(ComplexF64, size(tbt)), Val(C)),
) where {D, C}
    block = tbt.block
    block_i, block_j = tbt.block_ij
    is = tbt.axis[Block(block_i)] # global row indices
    js = tbt.axis[Block(block_j)] # global col indices
    MmtC = block.MmtC # contracted product of `Mm` and (complexified) `t`

    all(∈(1:D), components) || error(lazy"components must be in 1:$D; got $components")

    # NB: ↓ one more case of assuming no free parameters in 
    #     `δs = constant.(orbit(block.h_orbit))` (we don't materialize `δs` though)
    δs = orbit(block.h_orbit)
    v_conj = cispi.(dot.(Ref(-2 .* k), constant.(δs)))
    δ_mult_v_conj = similar(v_conj) # preallocate for reuse below
    for (idx, component) in enumerate(components)
        ∇H = ∇Hs[idx]
        # compute `δ[component]` multiplied by `v_conj` efficiently; equivalent to
        # `v_conj .* getindex.(constant.(δs), component)` but without repeated allocations
        unsafe_copyto!(δ_mult_v_conj, 1, v_conj, 1, length(v_conj))
        δ_mult_v_conj .*= getindex.(constant.(δs), component) # not yet including (-2πi) factor
        for (local_i, i) in enumerate(is)
            for (local_j, j) in enumerate(js)
                ∇Hᵢⱼ = @inbounds dot(δ_mult_v_conj, @view MmtC[:, local_i, local_j])
                ∇Hᵢⱼ *= 2im * π # include (-2πi) factor; now no minus sign (outside `dot`)
                isnothing(c) || (∇Hᵢⱼ *= c) # multiply by coefficient if provided
                ∇H[i, j] += ∇Hᵢⱼ
                i == j && continue # don't add diagonal elements twice
                ∇H[j, i] += tbt.hermiticity == ANTIHERMITIAN ? -conj(∇Hᵢⱼ) : conj(∇Hᵢⱼ)
            end
        end
    end

    return ∇Hs
end