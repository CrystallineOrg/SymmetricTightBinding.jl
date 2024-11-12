"""
Obtains directly the symmetry vectos for the bands computed in the MPB model `ms` for the space
group defined in `sg_num`. It fixs up the symmetry content at Γ and ω=0 and returns the symmetry
vectors and topoligies of the bands.
"""
function obtain_symmetry_vectors(ms::PyObject, sg_num::Int)
    brs = bandreps(sg_num) # elementary band representations
    lgs = littlegroups(sg_num) # little groups
    filter!(((klab, _),) -> klab ∈ klabels(brs), lgs) # restrict to k-points in `brs`
    map!(lg -> primitivize(lg, false), values(lgs)) # convert to primitive setting
    lgirsd = pick_lgirreps(lgs; timereversal=true) # small irreps associated with `lgs`

    symeigsd = Dict{String,Vector{Vector{ComplexF64}}}()
    for (klab, lg) in lgs
        kv = mp.Vector3(position(lg)()...)
        ms.solve_kpoint(kv)

        symeigsd[klab] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:ms.num_bands]
        for (i, gᵢ) in enumerate(lg)
            W = mp.Matrix(eachcol(rotation(gᵢ))...) # decompose gᵢ = {W|w}
            w = mp.Vector3(translation(gᵢ)...)
            symeigs = ms.compute_symmetries(W, w) # compute ⟨Eₙₖ|gᵢDₙₖ⟩ for all bands
            setindex!.(symeigsd[klab], symeigs, i) # update container of symmetry eigenvalues
        end
    end

    # --- fix singular photonic symmetry content at Γ, ω=0 --- 
    fixup_gamma_symmetry!(symeigsd, lgs)

    # --- analyze connectivity and topology of symmetry data ---
    summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)

    # --- convert to `SymmetryVector`s ---
    c_brs = calc_bandreps(sg_num, Val(3))
    symvecs = bandsum2symvec.(summaries, Ref(c_brs))
    topologies = getfield.(summaries, Ref(:topology))
    return symvecs, topologies
end

#= 
`t` -> dimension os the auxiliary modes to search
`brs` -> collection of the BRs of the SG
=#

"""
Finds all sets of bands in the SG that have dimension equal to `μᴸ`.

1. `μᴸ` -> dimension of the auxiliary modes to search
2. `brs` -> collection of the BRs of the SG
"""
function find_auxiliary_modes(μᴸ::Int, brs::Collection{<:NewBandRep})
    iszero(μᴸ) && return [Int[]]
    μs_brs = occupation.(brs)
    long_cand = find_all_admissible_expansions(
        brs, μs_brs, μᴸ, #= occupation =#
        Int[], Int[]) #= idxs =#

    return long_cand
end

"""
Computes the generalized inverse `Xᵍ` of `X`, computed from the Smith normal form.
"""
function generalized_inv(X::AbstractMatrix{<:Integer})
    F = smith(X)
    Λ = MPBUtils.diagm(F)
    Λg = zeros(Float64, size(Λ)[2], size(Λ)[1])
    for (n, λₙ) in enumerate(F.SNF)
        Λg[n, n] = iszero(λₙ) ? λₙ : inv(λₙ)
    end
    Xᵍ = F.Tinv * Λg * F.Sinv # generalized inverse

    return Xᵍ
end

#=
Different notation used explained here. First we define the notation for symmetry vectors
obtained from MPB vs the ones for the solutions:

    m                        =====> MPB
    nᴸ, nᵀ⁺ᴸ, nᵀ = nᵀ⁺ᴸ - nᴸ =====> solutions

Then, symmetry vectors can be explitted in several ways depending of if the irreps belong to
Γ or not and if the irreps belongs to higher frecuency bands or just ω=0:

    m = mᵧ + m₋ᵧ =====> Diffrentiate from Γ and not Γ
    n = nᵧ + n₋ᵧ =====> Diffrentiate from Γ and not Γ

    mᵧ = mᵧ⁼⁰ + mᵧꜛ⁰ =====> Diffrentiate from ω=0 and ω>0
    nᵧ = nᵧ⁼⁰ + nᵧꜛ⁰ =====> Diffrentiate from ω=0 and ω>0

with the notation clear, now we need to check if the solution we have obtained is physical
or not. In other words, we need to check two things:

1. Whether if our solution subduce properly the $O(3)$ representation at $\Gamma$ and zero 
frequency. This can be check easily using `PhotonicBandConnectivity.jl`. As estipulated 
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
Check if a certain solution `(nᵀ⁺ᴸ, nᴸ)` is *physical* given a symmetry vector `m`. A solution
is physical if it fulfills to checks:

    1. It subduces properly the O(3) representation at Γ and zero frequency.
    2. It doesn't make use of the higher frequency irreps to regularize the symmetry content 
    at zero frequency, and that instead uses the auxiliary modes `nᴸ` to achieve it.
"""
function is_integer_p_check(m::AbstractSymmetryVector,
    nᵀ⁺ᴸ::AbstractSymmetryVector{D},
    nᴸ::AbstractSymmetryVector{D},
    Q::Matrix{Int},
    Γidx::Int
) where {D}
    # convert everythin into vectors w/o occupation
    mᵧ = multiplicities(m)[Γidx]
    nᵀ⁺ᴸᵧ = multiplicities(nᵀ⁺ᴸ)[Γidx]
    nᴸᵧ = multiplicities(nᴸ)[Γidx]

    Q⁺ = generalized_inv(Q)
    nᵀᵧ = nᵀ⁺ᴸᵧ - nᴸᵧ # obtain the symmetry vector of the tranversal modes
    p = Q⁺ * (nᵀᵧ - mᵧ) # compute the vector p

    # finally check if the vector p is an integer vector and if all the irreps with 
    # negative multiplicite are present on the longitudinal modes
    p_int = round.(Int, p)
    p_int ≈ p || error("unexpectedly found non-integer p - unhandled")
    return p_int
end

"""
Obtains a possible TETB model `nᵀ⁺ᴸ` for the auxiliary modes provided `idxsᴸs`.
"""
function find_apolar_modes(
    m::AbstractSymmetryVector{D},
    idxsᴸs::Vector{Vector{Int64}},
    brs::Collection{NewBandRep{D}}) where {D}

    μs_brs = occupation.(brs)
    idxs = eachindex(first(brs))

    # compute the fixed part and the free part of the physical ω=0 irreps at Γ
    Γidx = something(findfirst(==("Γ"), klabels(m)))
    lgirs = irreps(m)[Γidx]

    nfixed, Q = physical_zero_frequency_gamma_irreps(
        lgirs;
        supergroup_constraints=true,
        force_fixed=true,
        lattice_reduce=true)

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
        #   @Γ : nᴸ[i] ≥ -nfixed ==> nᵀ⁺ᴸ[i] ≥ (m - n_fixed)[i]
        # We can fold these two sets of constraints into one, via the following 
        # manipulations:
        constraints = m + nᴸ # now the constraints are wrong at Γ; proceed to correct this
        constraints.multsv[Γidx] -= nfixed + multiplicities(nᴸ)[Γidx] # now: fixed

        idxsᵀ⁺ᴸs = find_all_admissible_expansions(brs, μs_brs, μᵀ⁺ᴸ, Vector(constraints), idxs)

        if !isempty(idxsᵀ⁺ᴸs)
            ps = map(idxsᵀ⁺ᴸs) do idxsᵀ⁺ᴸ
                nᵀ⁺ᴸ = SymmetryVector(sum(brs[idxsᵀ⁺ᴸ]))
                is_integer_p_check(m, nᵀ⁺ᴸ, nᴸ, Q, Γidx)
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
made. This parameter is settle by default to `true`
"""
function find_bandrep_decompositions(
    m::AbstractSymmetryVector{D},
    brs::Collection{NewBandRep{D}};
    μᴸ_min::Integer=0,
    μᴸ_max::Integer=μᴸ_min + 2 * occupation(m),
    connected_to_zero_frequency::Bool=true
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

# we do not include the (usually redundant) exponential (k-dependet) phases below
"""
Induce a representation for the generators of the SG from a representation of the site-symmetry 
group of a particular maximal WP.
"""
function sgrep_induced_by_siteir_generators(br::NewBandRep{D}) where {D}
    siteir = br.siteir
    siteir_dim = irdim(siteir)
    siteg = group(siteir)
    wps = orbit(siteg)
    mult = length(wps)
    gens = generators(num(siteg), SpaceGroup{D})

    ρs = [BlockArray{ComplexF64}(
        zeros(ComplexF64, siteir_dim * mult, siteir_dim * mult),
        fill(siteir_dim, mult), fill(siteir_dim, mult)) for _ in eachindex(gens)]
    for (n, g) in enumerate(gens)
        ρ = ρs[n]
        for (α, (gₐ, qₐ)) in enumerate(zip(cosets(siteg), wps))
            check = false
            for (β, (gᵦ, qᵦ)) in enumerate(zip(cosets(siteg), wps))
                tᵦₐ = constant(g * parent(qₐ) - parent(qᵦ)) # ignore free parts of the WP
                # compute h = gᵦ⁻¹ tᵦₐ⁻¹ g gₐ
                h = compose(compose(compose(inv(gᵦ), SymOperation(-tᵦₐ), false), g, false), gₐ, false)
                idx_h = findfirst(==(h), siteg)
                if !isnothing(idx_h) # h ∈ siteg and qₐ and qᵦ are connected by g
                    ρ[Block(β, α)] .= siteir.matrices[idx_h]
                    check = true
                    break
                end
            end
            check || error("failed to find any nonzero block")
        end
    end

    return gens .=> ρs
end

"""
Direct sum of matrices
"""
function ⊕(As::AbstractMatrix...)
    return cat(As..., dims=Val((1, 2)))
end


function sgrep_induced_by_siteir_generators(brs::CompositeBandRep{D}) where {D}
    gens = generators(num(brs), SpaceGroup{D})
    ρs = [Matrix{Complex}(undef, 0, 0) for _ in eachindex(gens)]

    for (idxc, c) in enumerate(brs.coefs)
        sgrep = sgrep_induced_by_siteir_generators(brs.brs[idxc])
        for idxg in eachindex(gens)
            iszero(c) && continue
            ρ_idxg = sgrep[idxg][2]
            for _ in 1:Int(c)
                ρs[idxg] = ρs[idxg] ⊕ ρ_idxg
            end
        end
    end
    return gens .=> ρs
end

"""
Compute the symmetry related hoppings relative vectors from the WP of `br1` to the WP of `br2`
displaced a set of primitive lattice vectors `Rs`.
"""
function obtain_symmetry_related_hoppings(
    Rs::AbstractVector{V}, # must be specified in the primitive basis
    br1::NewBandRep{D},
    br2::NewBandRep{D}
) where {V<:Union{AbstractVector{<:Integer},RVec{D}}} where {D}

    sgnum = num(br1)
    num(br2) == sgnum || error("both band representations must be in the same space group")
    # we only want to include the wyckoff positions in the primitive cell - but the default
    # listings from `spacegroup` include operations that are "centering translations";
    # fortunately, the orbit returned for a `NewBandRep` do not include these redundant
    # operations - but is still specified in a conventional basis. So, below, we remove
    # redundant operations from the space group, and also change both the operations and the
    # positions from a conventional to a primitive basis
    cntr = centering(sgnum, D)
    ops = primitivize(spacegroup(sgnum, D))
    wps1 = primitivize.(orbit(group(br1)), cntr)
    wps2 = primitivize.(orbit(group(br2)), cntr)

    # we are going to create a dictionary of dictionaries. The first set of keys will indicate
    # the representative of the hopping distance and inside it, each dictionary will be 
    # associated with a point in the "orbit" of that hopping distance. Addtionally, if you 
    # look at the value of and specific point inside the "orbit" of one representative you 
    # will obtain a set of tuples that will encode the points of the WPs involved in the 
    # hopping and the displacement vector R. (q_i, w_j , R) => q_i -> w_j + R

    orbsdd = Dict{RVec{D},Dict{RVec{D},Vector{Tuple{RVec{D},RVec{D},RVec{D}}}}}()
    for R in Rs
        R = RVec(R) # change the type of R to be type consistent
        for (qₐ, qᵦ) in Iterators.product(wps1, wps2)
            qₐ = parent(qₐ) # to work with the RVec type directly
            qᵦ = parent(qᵦ) # same as above
            δ = qᵦ + R - qₐ
            if !in(δ, keys(orbsdd)) && !in(δ, Iterators.flatten(keys.(values(orbsdd))))
                push!(orbsdd, δ => Dict(δ => [(qₐ, qᵦ, R)]))
                for g in ops
                    Ρ = SymOperation(rotation(g)) # type consitentency for the rotation 
                    # we want to keep the WPs inside the unit cell and put all the posible
                    # translations into the RVec 'R', so we can keep the notation such that
                    # q_i -> w_j + R. For that reason we make the following changes:
                    qₐ′ = RVec(Crystalline.reduce_translation_to_unitrange(Crystalline.constant(Ρ * qₐ)))
                    qᵦ′ = RVec(Crystalline.reduce_translation_to_unitrange(Crystalline.constant(Ρ * qᵦ)))
                    dₐ = (Ρ * qₐ) - qₐ′
                    dᵦ = (Ρ * qᵦ) - qᵦ′
                    R′ = (Ρ * R) + dᵦ - dₐ
                    δ′ = Ρ * δ
                    if !in(δ′, keys(orbsdd)) && !in(δ′, keys(orbsdd[δ]))
                        push!(orbsdd[δ], δ′ => [(qₐ′, qᵦ′, R′)])
                    elseif in(δ′, keys(orbsdd[δ])) && !in((qₐ′, qᵦ′, R′), orbsdd[δ][δ′])
                        push!(orbsdd[δ][δ′], (qₐ′, qᵦ′, R′))
                    else
                        continue
                    end
                end
            else
                continue
            end
        end
    end
    return orbsdd
end