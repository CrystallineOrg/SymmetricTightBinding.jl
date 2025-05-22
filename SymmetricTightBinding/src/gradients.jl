struct TightBindingModelGradient{D}
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
gradient_wrt_hopping(tbm::TightBindingModel) = TightBindingModelGradient(tbm)
gradient_wrt_hopping(ptbm::ParameterizedTightBindingModel) = gradient_wrt_hopping(ptbm.tbm)

function (tbmg::TightBindingModelGradient{D})(k::ReciprocalPointLike{D}) where D
    map(tbt->tbt(k), tbmg.tbm.terms)
end

# evaluate just the `i`th component of the coefficient gradient
function (tbmg::TightBindingModelGradient{D})(k::ReciprocalPointLike{D}, i::Int) where D
    tbmg.tbm[i](k)
end

"""
    energy_gradient_wrt_hopping(
        ptbm::ParameterizedTightBindingModel{D},
        k::ReciprocalPointLike{D};
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
    k::ReciprocalPointLike{D};
    degen_rtol::Float64 = 1e-12,
    degen_atol::Float64 = 1e-12
) where D
    Nᶜ = length(ptbm.tbm) # number of hopping terms
    Nᵇ = ptbm.tbm.N       # number of bands

    Es, us = solve(ptbm, k; bloch_phase=Val(false)) # "unperturbed" energies and eigenstates
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
