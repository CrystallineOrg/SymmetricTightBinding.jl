"""
Erase the content of a certain HSP from the input object. This HSP is settle by default to Γ.
"""
function prune_klab_irreps!(brs::Collection{<:NewBandRep}, klab::String="Γ")
    prune_kidx = findfirst(==(klab), klabels(brs))
    isnothing(prune_kidx) && error(lazy"could not find $klab among included k-points")

    lgirsv = irreps(first(brs))
    foreach(brs) do br
        deleteat!(multiplicities(br.n), prune_kidx)
        @assert irreps(br.n) === lgirsv
    end
    deleteat!(lgirsv, prune_kidx)

    return brs
end

"""
Copies the symmetry vectors of the bandreps inside the collection so they are not overwritten
"""
function _symmetry_vector_shallow_copy(brs::Collection{<:NewBandRep})
    lgirsv′ = copy(irreps(first(brs)))
    brs′ = Collection(map(brs) do br
        NewBandRep(
            br.siteir,
            SymmetryVector(lgirsv′, copy(multiplicities(br)), occupation(br)),
            br.spinful,
            br.timereversal)
    end)
    return brs′
end
"""
Erase the content of a certain HSP from the input object. This HSP is settle by default to Γ.
"""
function prune_klab_irreps(brs::Collection{<:NewBandRep}, klab::String="Γ")
    return prune_klab_irreps!(_symmetry_vector_shallow_copy(brs), klab)
end

function prune_klab_irreps!(v::SymmetryVector, klab::String="Γ")
    prune_iridxs = findall(lgirs -> klabel(first(lgirs)) == klab, irreps(v))
    isempty(prune_iridxs) && error(lazy"could not find $klab among included irreps")
    deleteat!(multiplicities(v), prune_iridxs)
    deleteat!(irreps(v), prune_iridxs)
    return v
end
function prune_klab_irreps!(v::AbstractSymmetryVector, klab::String="Γ")
    prune_klab_irreps!(SymmetryVector(v), klab)
end

function prune_klab_irreps(v::SymmetryVector, klab::String="Γ")
    lgirsv´ = copy(irreps(v))
    multsv´ = copy(multiplicities(v))
    v´ = SymmetryVector(lgirsv´, multsv´, occupation(v))
    return prune_klab_irreps!(v´, klab)
end
function prune_klab_irreps(v::AbstractSymmetryVector, klab::String="Γ")
    prune_klab_irreps(SymmetryVector(v), klab)
end

"""
Picks the content of a certain HSP from the input object, erasing all the content for other 
HSPs. This HSP is settle by default to Γ.
"""
function pick_klab_irreps!(brs::Collection{<:NewBandRep}, klab::String="Γ")
    prune_kidx = findall(!=(klab), klabels(brs))
    isnothing(prune_kidx) && error(lazy"could not find $klab among included k-points")

    lgirsv = irreps(first(brs))
    foreach(brs) do br
        deleteat!(multiplicities(br.n), prune_kidx)
        @assert irreps(br.n) === lgirsv
    end
    deleteat!(lgirsv, prune_kidx)

    return brs
end

function pick_klab_irreps(brs::Collection{<:NewBandRep}, klab::String="Γ")
    return pick_klab_irreps!(_symmetry_vector_shallow_copy(brs), klab)
end


function pick_klab_irreps!(v::SymmetryVector, klab::String="Γ")
    prune_iridxs = findall(lgirs -> klabel(first(lgirs)) != klab, irreps(v))
    isempty(prune_iridxs) && error(lazy"could not find $klab among included irreps")
    deleteat!(multiplicities(v), prune_iridxs)
    deleteat!(irreps(v), prune_iridxs)
    return v
end
function pick_klab_irreps!(v::AbstractSymmetryVector, klab::String="Γ")
    pick_klab_irreps!(SymmetryVector(v), klab)
end


function pick_klab_irreps(v::SymmetryVector, klab::String="Γ")
    lgirsv´ = copy(irreps(v))
    multsv´ = copy(multiplicities(v))
    v´ = SymmetryVector(lgirsv´, multsv´, occupation(v))
    return pick_klab_irreps!(v´, klab)
end
function pick_klab_irreps(v::AbstractSymmetryVector, klab::String="Γ")
    pick_klab_irreps(SymmetryVector(v), klab)
end

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

2. Whether if our solution doesn't make use of the higher frequency irreps present in 
$m_\Gamma^{>0}$ to regularize the symmetry content at zero frequency, and that instead those 
negative multiplicities in the irreps are cancel out by the longitudinal modes $n^L$. This 
is checked in the following way:

    - If $n_{\Gamma,i}^{T+L} > m^{>0}_\Gamma \quad \forall i$, all higher frequency irreps 
    are present on the model so the model perfectly represent the system bans.
    - If $n_{\Gamma,i}^{T+L} < m^{>0}_\Gamma$ for some $i$ that higher frequency irrep is 
    lacking in the model so the solution do not proper represent the higher frequency bands 
    so it must be discarded.

    Those conditions can be easily checked by the following condition:

    $$ n_{\Gamma,i}^{T+L} - m^{>0}_\Gamma = abs(n_{\Gamma,i}^{T+L} - m^{>0}_\Gamma) $$ 

=#

"""
Check if a certain solution `(nᵀ⁺ᴸ, nᴸ)` is *physical* given a symmetry vector `m`. A solution
is physical if it fulfills to checks:

    1. It subduces properly the O(3) representation at Γ and zero frequency.
    2. It doesn't make use of the higher frequency irreps to regularize the symmetry content 
    at zero frequency, and that instead uses the auxiliary modes `nᴸ` to cancel them.
"""
function is_integer_p_check(m::AbstractSymmetryVector,
        nᵀ⁺ᴸ::AbstractSymmetryVector,
        nᴸ::AbstractSymmetryVector,
        Q::Matrix{Int},
        Γidx::Int
    )
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
    p_int≈p || error("unexpectedly found non-integer p - unhandled")
    return p_int
end

"""
Obtains a possible TETB model `nᵀ⁺ᴸ` for the auxiliary modes provided `idxsᴸs`. Additionally,
it checks if the solution provided is physical or not.
"""
function find_apolar_modes(
    m::AbstractSymmetryVector,
    idxsᴸs::Vector{Vector{Int64}},
    brs::Collection{<:NewBandRep})

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
        #   @Γ : nᵀ⁺ᴸ[i] ≥ (m - n_fixed)[i]
        #   @-Γ: nᵀ⁺ᴸ[i] == (m + nᴸ)[i]   (and we translate this to nᵀ⁺ᴸ[i] ≥ (m + nᴸ)[i]
        #                                  cf. non-negativity)
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

            candidates = TightBindingCandidateSet(idxsᴸ, idxsᵀ⁺ᴸs, ps, brs)
            push!(candidatesv, candidates)
        end
    end
    return candidatesv
end


function find_bandrep_decompositions(
        m::AbstractSymmetryVector{D},
        brs::Collection{NewBandRep{D}};
        μᴸ_min :: Integer = 0,
        μᴸ_max :: Integer = μᴸ_min + 2*occupation(m),
        connected_to_zero_frequency :: Bool = true
    ) where D

    connected_to_zero_frequency || error("not implemented yet") # TODO

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