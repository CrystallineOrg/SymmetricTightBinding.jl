# Up-front caching of the Hamiltonian terms hᵢ(k) over a fixed set of k-points.
#
# The Hamiltonian is linear in the hopping coefficients, H(k) = ∑ᵢ csᵢhᵢ(k), with the term
# matrices hᵢ(k) independent of the coefficients `cs`. Workloads that repeatedly evaluate a
# model at the same k-points but varying coefficients — fitting, above all — therefore
# needlessly recompute every hᵢ(k) (Bloch phases and all) on every evaluation, twice over if
# they also need coefficient gradients. `TightBindingCache` removes this redundancy by
# tabulating hᵢ(k) once, up front, reducing Hamiltonian assembly to an accumulation of
# cached matrices and the Feynman–Hellmann coefficient gradient to inner products against
# them; both operations write into preallocated work arrays, avoiding allocation entirely.

"""
    TightBindingCache(tbm::TightBindingModel{D}, ks::AbstractVector{<:ReciprocalPointLike{D}})

A cache of the coefficient-independent Hamiltonian term matrices ``h_i(\\mathbf{k})`` of
`tbm` (i.e., ``H(\\mathbf{k}) = \\sum_i c_i h_i(\\mathbf{k})``), tabulated up front over the
fixed **k**-points `ks`, plus preallocated work arrays for Hamiltonian assembly and for
Feynman–Hellmann coefficient gradients. Intended for workloads — chiefly, fitting — that
evaluate the same model over the same **k**-points many times with varying coefficients.

The memory footprint is dominated by the term table: `length(tbm) * length(ks)` complex
matrices of size `tbm.N × tbm.N`.

## Usage

With `cache = TightBindingCache(tbm, ks)` and coefficients `cs`:
- `cache(cs, κ)`: assemble the Hamiltonian at `ks[κ]`, i.e., the cached analogue of
  `tbm(cs)(ks[κ])`. The result is written into (and aliases) the internal work array
  `cache.W`; it may be freely mutated, e.g., by `eigen!`/`eigvals!`, and is overwritten by
  the next assembly call.
- `energy_gradient_wrt_hopping(cache, κ, (Es, us))`: the cached analogue of
  `energy_gradient_wrt_hopping(tbm(cs), ks[κ], (Es, us))`. The result aliases the internal
  buffer `cache.∇ᶜEs` and is overwritten by the next gradient call.

Because evaluations share the internal work arrays, a cache must not be used from multiple
threads concurrently (create one cache per thread instead; only the work arrays are then
duplicated meaningfully).
"""
struct TightBindingCache{D, S, T, K <: AbstractVector{<:ReciprocalPointLike{D}}}
    tbm  :: TightBindingModel{D, S}
    ks   :: K                                   # the fixed k-points of the cache
    hs   :: Vector{Vector{Matrix{ComplexF64}}}  # hs[κ][i] = hᵢ(ks[κ])
    W    :: Matrix{ComplexF64}                  # assembly work array (N×N); see docstring
    ∇ᶜEs :: Matrix{T}                           # Feynman–Hellmann gradient buffer (Nᶜ×Nᵇ)
end

function TightBindingCache(
    tbm::TightBindingModel{D, S},
    ks::AbstractVector{<:ReciprocalPointLike{D}},
) where {D, S}
    Nᶜ, N = length(tbm), tbm.N
    hs = [[tbm[i](k) for i in 1:Nᶜ] for k in ks]
    W = Matrix{ComplexF64}(undef, N, N)
    T = S === HERMITIAN ? Float64 : ComplexF64 # cf. `energy_gradient_wrt_hopping`
    ∇ᶜEs = Matrix{T}(undef, Nᶜ, N)
    return TightBindingCache{D, S, T, typeof(ks)}(tbm, ks, hs, W, ∇ᶜEs)
end

# Hamiltonian assembly at the `κ`th cached k-point: H(ks[κ]) = ∑ᵢ csᵢhᵢ(ks[κ]), accumulated
# into the work array `cache.W` (cf. `TightBindingCache`'s docstring on aliasing/mutation)
function (cache::TightBindingCache{D, S})(
    cs::AbstractVector{<:Real},
    κ::Integer,
) where {D, S}
    length(cs) ≠ length(cache.tbm) && _throw_term_coef_length_mismatch(cache.tbm.terms, cs)
    W = fill!(cache.W, zero(ComplexF64))
    for (c, h) in zip(cs, cache.hs[κ])
        W .+= c .* h
    end
    return S == HERMITIAN ? Hermitian(W) : W
end

"""
    energy_gradient_wrt_hopping(
        cache::TightBindingCache{D},
        κ::Integer,
        (Es, us);
        degen_rtol::Float64 = 1e-12,
        degen_atol::Float64 = 1e-12
    ) where D

Cached variant of `energy_gradient_wrt_hopping(ptbm, k, (Es, us))`, evaluated at the `κ`th
cached k-point of `cache` using its tabulated term matrices hᵢ(ks[κ]): the eigensolution
`(Es, us)` must correspond to `ks[κ]` (e.g., obtained by diagonalizing `cache(cs, κ)`).

The returned column views alias the cache's internal buffer `cache.∇ᶜEs` and are overwritten
by the next gradient call.
"""
function energy_gradient_wrt_hopping(
    cache::TightBindingCache{D, S},
    κ::Integer,
    (Es, us);
    degen_rtol::Float64 = 1e-12,
    degen_atol::Float64 = 1e-12,
) where {D, S}
    if S === NONHERMITIAN
        error("energy gradient with respect to hopping is not currently implemented for \
               NONHERMITIAN models")
    end
    bands = _degenerate_band_groups(Es, degen_rtol, degen_atol)

    # apply the Feynman–Hellmann theorem along each hopping coefficient `i`
    ∇ᶜEs = cache.∇ᶜEs
    for (i, ∂H) in enumerate(cache.hs[κ])
        _feynman_hellmann_derivatives!((@view ∇ᶜEs[i, :]), ∂H, us, bands, Val(S))
    end

    return eachcol(∇ᶜEs)
end

# (kept here, rather than in show.jl, since `caching.jl` is included after `show.jl`)
function Base.show(io::IO, ::MIME"text/plain", cache::TightBindingCache{D}) where {D}
    Nᵏ = length(cache.ks)
    print(io, "TightBindingCache{", D, "} over ", Nᵏ, " k-point", Nᵏ == 1 ? "" : "s", ":\n ")
    summary(io, cache.tbm)
end
