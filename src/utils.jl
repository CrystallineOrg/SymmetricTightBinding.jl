using Crystalline: translation

"""
    obtain_symmetry_vectors(ms::PyObject, sgnum::Int) --> Tuple{Vector{SymmetryVector{D}}, Vector{TopologyKind}}

Obtains directly the symmetry vector for the bands computed in the MPB model `ms` for the space
group defined in `sgnum`. It fixes up the symmetry content at Γ and ω=0 and returns the symmetry
vectors and topologies of the bands.
"""
function obtain_symmetry_vectors(ms::PyObject, sgnum::Int)
    brs = bandreps(sgnum) # elementary band representations
    lgs = littlegroups(sgnum) # little groups
    filter!(((klab, _),) -> klab ∈ klabels(brs), lgs) # restrict to k-points in `brs`
    map!(lg -> primitivize(lg, false), values(lgs)) # convert to primitive setting
    lgirsd = pick_lgirreps(lgs; timereversal) # small irreps associated with `lgs`

    symeigsd = Dict{String, Vector{Vector{ComplexF64}}}()
    for (klab, lg) in lgs
        kv = mp.Vector3(position(lg)()...)
        ms.solve_kpoint(kv)

        symeigsd[klab] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:ms.num_bands]
        for (i, gᵢ) in enumerate(lg)
            W = mp.Matrix(eachcol(rotation(gᵢ))...) # decompose gᵢ = {W|w}
            w = mp.Vector3(translation(gᵢ)...)
            symeigs = ms.compute_symmetries(W, w) # compute ⟨Eₙₖ|gᵢ Dₙₖ⟩ for all bands
            setindex!.(symeigsd[klab], symeigs, i) # update container of symmetry eigenvalues
        end
    end

    # --- fix singular photonic symmetry content at Γ, ω=0 --- 
    fixup_gamma_symmetry!(symeigsd, lgs)

    # --- analyze connectivity and topology of symmetry data ---
    summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)

    # --- convert to `SymmetryVector`s ---
    c_brs = calc_bandreps(sgnum, Val(3))
    symvecs = bandsum2symvec.(summaries, Ref(c_brs))
    topologies = getfield.(summaries, Ref(:topology))
    return symvecs, topologies
end

#= 
`t` -> dimension os the auxiliary modes to search
`brs` -> collection of the BRs of the SG
=#

"""
    find_auxiliary_modes(μᴸ::Int, brs::Collection{<:NewBandRep}) -> Vector{Vector{Int}}

Finds all sets of bands in the SG that have dimension equal to `μᴸ`.

1. `μᴸ` -> dimension of the auxiliary modes to search
2. `brs` -> collection of the BRs of the SG
"""
function find_auxiliary_modes(μᴸ::Int, brs::Collection{<:NewBandRep})
    iszero(μᴸ) && return [Int[]]
    μs_brs = occupation.(brs)
    long_cand = find_all_admissible_expansions(
        brs,
        μs_brs,
        μᴸ, #= occupation =#
        Int[],
        Int[],
    ) #= idxs =#

    return long_cand
end

"""
    generalized_inv(X::AbstractMatrix{<:Integer}) -> AbstractMatrix{Float64}

Computes the generalized inverse `Xᵍ` of `X`, computed from the Smith normal form.
"""
function generalized_inv(X::AbstractMatrix{<:Integer})
    F = smith(X) # compute the Smith normal form X = SΛT
    Λ = MPBUtils.diagm(F)
    Λg = zeros(Float64, size(Λ)[2], size(Λ)[1])
    for (n, λₙ) in enumerate(F.SNF)
        Λg[n, n] = iszero(λₙ) ? λₙ : inv(λₙ) # inverse of Λ considering the zero values
    end
    Xᵍ = F.Tinv * Λg * F.Sinv # generalized inverse

    return Xᵍ
end

#=
Different notation used explained here. First we define the notation for symmetry vectors
obtained from MPB vs the ones for the solutions:

    m                        =====> MPB
    nᴸ, nᵀ⁺ᴸ, nᵀ = nᵀ⁺ᴸ - nᴸ =====> solutions

Then, symmetry vectors can be split in several ways depending of if the irreps belong to
Γ or not and if the irreps belongs to higher frequency bands or just ω=0:

    m = mᵧ + m₋ᵧ =====> Differentiate from Γ and not Γ
    n = nᵧ + n₋ᵧ =====> Differentiate from Γ and not Γ

    mᵧ = mᵧ⁼⁰ + mᵧꜛ⁰ =====> Differentiate from ω=0 and ω>0
    nᵧ = nᵧ⁼⁰ + nᵧꜛ⁰ =====> Differentiate from ω=0 and ω>0

with the notation clear, now we need to check if the solution we have obtained is physical
or not. In other words, we need to check two things:

1. Whether if our solution subduces properly the $O(3)$ representation at $\Gamma$ and zero 
frequency. This can be check easily using `PhotonicBandConnectivity.jl`. As stipulated 
before in [Problem 2](#problem-2), this is fulfilled if $\mathbf{p}\in\mathbb{Z}$.

2. Whether our solution doesn't make use of the higher frequency irreps present in 
$m_\Gamma^{>0}$ to regularize the symmetry content at zero frequency, and that instead 
those negative multiplicities in the irreps are cancelled out by the longitudinal modes $n^L$. 
We ensure this by the following check:

    Define the candidate-solution's zero-frequency content at $\Gamma$ by:

    $$n_\Gamma^{T,=0} = n_{\Gamma}^{T} - m_{\Gamma}^{>0} = n_{\Gamma}^{T+L} - n_{\Gamma}^L 
    - m_{\Gamma}^{>0} = m_{\Gamma}^{=0} + Q\mathbf{p}.$$

    Consider the following two cases:
    - If $n_{\Gamma,i}^{T,=0} < 0$ for some $i$, then $n_{\Gamma,i}^L \geq |n_{\Gamma,i}^{T,=0}|$ 
        for that $i$; equivalently, in this case $n_{\Gamma,i}^L \geq -n_{\Gamma,i}^{T,=0}$.
    - Conversely, if  $n_{\Gamma,i}^{T,=0} ≥ 0$ for some $i$, we still have $n_{\Gamma,i}^L ≥ 0$
         and consequently also $n_{\Gamma,i}^L ≥ -n_{\Gamma,i}^{T,=0}$.

    Thus, regardless of the sign of $n_{\Gamma,i}^{T,=0}$, we may require that:

    $$ n_{\Gamma}^L \geq -n_\Gamma^{T,=0}$$

=#

"""
            is_integer_p_check(m::AbstractSymmetryVector,
                nᵀ⁺ᴸ::AbstractSymmetryVector{D},
                nᴸ::AbstractSymmetryVector{D},
                Q::Matrix{Int},
                Γ_idx::Int) where {D}

Check if a certain solution `(nᵀ⁺ᴸ, nᴸ)` subduces properly the O(3) representation at Γ and 
zero frequency. This is fulfilled if `p` is an integer vector. Check `devdocs.md`.
"""
function is_integer_p_check(
    m::AbstractSymmetryVector,
    nᵀ⁺ᴸ::AbstractSymmetryVector{D},
    nᴸ::AbstractSymmetryVector{D},
    Q::Matrix{Int},
    Γ_idx::Int,
) where {D}
    # convert everything into vectors w/o occupation and take the content at Γ
    mᵧ = multiplicities(m)[Γ_idx]
    nᵀ⁺ᴸᵧ = multiplicities(nᵀ⁺ᴸ)[Γ_idx]
    nᴸᵧ = multiplicities(nᴸ)[Γ_idx]

    Q⁺ = generalized_inv(Q)
    nᵀᵧ = nᵀ⁺ᴸᵧ - nᴸᵧ # obtain the symmetry vector of the transversal modes
    p = Q⁺ * (nᵀᵧ - mᵧ) # compute the vector p
    # this way we assure that the auxiliary modes are used to regularize the symmetry content
    # at zero frequency

    # finally check if the vector p is an integer vector
    p_int = round.(Int, p)
    p_int ≈ p || error("unexpectedly found non-integer p - unhandled")
    return p_int
end

"""
    find_apolar_modes(m::AbstractSymmetryVector{D},
                      idxsᴸs::Vector{Vector{Int}}, 
                      brs::Collection{NewBandRep{D}}) -> Vector{TightBindingCandidateSet}

Obtains a possible TETB model `nᵀ⁺ᴸ` for the auxiliary modes provided `idxsᴸs`.
    
It checks its **physicality** by ensuring that the solution subduces properly the O(3) 
representation at Γ and zero frequency and that the higher frequency irreps present in `m` 
are not used to regularize the symmetry content at zero frequency.
"""
function find_apolar_modes(
    m::AbstractSymmetryVector{D},
    idxsᴸs::Vector{Vector{Int64}},
    brs::Collection{NewBandRep{D}},
) where {D}
    μs_brs = occupation.(brs)
    idxs = eachindex(first(brs))

    # compute the fixed part `n_fixed` and the free part `Q` of the physical ω=0 irreps at Γ
    Γ_idx = something(findfirst(==("Γ"), klabels(m)))
    lgirs = irreps(m)[Γ_idx]

    n_fixed, Q = physical_zero_frequency_gamma_irreps(
        lgirs;
        supergroup_constraints = true,
        force_fixed = true,
        lattice_reduce = true,
    )

    candidatesv = TightBindingCandidateSet[]
    for idxsᴸ in idxsᴸs
        nᴸ = if isempty(idxsᴸ)
            zero(first(brs))
        else
            SymmetryVector(sum(brs[idxsᴸ]))
        end
        μᵀ⁺ᴸ = occupation(m) + occupation(nᴸ)

        # We want to enforce two constraints, one at Γ, one at "not-Γ" ≡ -Γ:
        #   @-Γ: nᵀ⁺ᴸ[i] == (m + nᴸ)[i]   (and we translate this to nᵀ⁺ᴸ[i] ≥ (m + nᴸ)[i]
        #                                  cf. non-negativity)
        #   @Γ : nᴸ[i] ≥ -n_fixed ==> nᵀ⁺ᴸ[i] ≥ (m - n_fixed)[i]
        # We can fold these two sets of constraints into one, via the following 
        # manipulations:
        constraints = m + nᴸ # now the constraints are wrong at Γ; proceed to correct this
        constraints.multsv[Γ_idx] -= n_fixed + multiplicities(nᴸ)[Γ_idx] # now: fixed

        idxsᵀ⁺ᴸs =
            find_all_admissible_expansions(brs, μs_brs, μᵀ⁺ᴸ, Vector(constraints), idxs)

        if !isempty(idxsᵀ⁺ᴸs)
            ps = map(idxsᵀ⁺ᴸs) do idxsᵀ⁺ᴸ
                nᵀ⁺ᴸ = SymmetryVector(sum(brs[idxsᵀ⁺ᴸ]))
                is_integer_p_check(m, nᵀ⁺ᴸ, nᴸ, Q, Γ_idx)
            end
            longitudinal = CompositeBandRep_from_indices(idxsᴸ, brs)
            apolarv = CompositeBandRep_from_indices.(idxsᵀ⁺ᴸs, Ref(brs))
            candidates = TightBindingCandidateSet(longitudinal, apolarv, ps)
            push!(candidatesv, candidates)
        end
    end
    return candidatesv
end

"""
Obtain a bandrep decomposition for the symmetry vector of the bands provided `m` with a minimal
number of auxiliary bands in the interval `[μᴸ_min,μᴸ_max]`.

If the photonic bands are connected to zero frequency corrections to the singularity at Γ are
made. This parameter is set by default to `true`.
"""
function find_bandrep_decompositions(
    m::AbstractSymmetryVector{D},
    brs::Collection{NewBandRep{D}};
    μᴸ_min::Integer = 0,
    μᴸ_max::Integer = μᴸ_min + 2 * occupation(m),
    connected_to_zero_frequency::Bool = true,
) where {D}
    if D < 3 || !connected_to_zero_frequency
        μ = occupation(m)
        μs_brs = occupation.(brs)
        idxs_k = eachindex(first(brs))

        # we don't need any longitudinal modes so we can directly find the expansions using
        # `m` as a constraint
        idxs_sol = find_all_admissible_expansions(brs, μs_brs, μ, Vector(m), idxs_k)

        longitudinal = Crystalline.CompositeBandRep_from_indices(Int[], brs)
        apolar = Crystalline.CompositeBandRep_from_indices.(idxs_sol, Ref(brs))

        !isempty(apolar) || error("Check the symmetry vector and space group used")

        return TightBindingCandidateSet(longitudinal, apolar, [Float64[]])

    else
        μᴸ = μᴸ_min - 1
        while μᴸ < μᴸ_max
            μᴸ += 1
            idxsᴸs = find_auxiliary_modes(μᴸ, brs)
            (isempty(idxsᴸs) && μᴸ ≠ 0) && continue
            # compute all possible decomposition of m into valid combinations of nᴸ and nᵀ⁺ᴸ
            candidatesv = find_apolar_modes(m, idxsᴸs, brs)
            isempty(candidatesv) || return candidatesv
        end
        error("""failed to find possible auxiliary-apolar decompositions for provided \
                symmetry vector in search range for auxiliary modes; increasing kwarg \
                `μᴸ_max` may help, if a decomposition exists""")
    end
end
