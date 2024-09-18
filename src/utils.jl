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
function prune_klab_irreps(brs::Collection{<:NewBandRep}, klab::String="Γ")
    return prune_klab_irreps!(_symmetry_vector_shallow_copy(brs), klab)
end

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

function find_auxiliary_modes(t::Int, d::Vector{Int64}, brs::Collection{<:NewBandRep})
    long_cand = find_all_admissible_expansions(
        brs, d, t, #= occupation =#
        Int[], Int[]) #= idxs =#

    return long_cand
end

# Computes the generalized inverse `Xᵍ` of `X`, computed from the Smith normal form.
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

We can obtain nᵀᵧ⁼⁰ from mᵧ⁼⁰ by:

    nᵀᵧ⁼⁰ = mᵧ⁼⁰ + Q*p

Now if p ∈ Ζ, the solution will be physical otherwise not. Additionally we must check if all
the irreps for ω>0 are reproduced so then:

    nᵀᵧꜛ⁰ - mᵧꜛ⁰ = nᵀᵧ - nᵀᵧ⁼⁰ - mᵧꜛ⁰ == 0

=#

function physical(mᵧ::AbstractSymmetryVector,
    nᵀ⁺ᴸᵧ::AbstractSymmetryVector,
    nᴸᵧ::AbstractSymmetryVector,
    nfree::Vector{Int},
    Q::Matrix{Int})
    # convert everythin into vectors w/o occupation
    mᵧ = Vector(mᵧ)[1:end-1]
    nᵀ⁺ᴸᵧ = Vector(nᵀ⁺ᴸᵧ)[1:end-1]
    nᴸᵧ = Vector(nᴸᵧ)[1:end-1]

    Q⁺ = generalized_inv(Q)
    nᵀᵧ = nᵀ⁺ᴸᵧ - nᴸᵧ # obtain the symmetry vector of the tranversal modes
    mᵧꜛ⁰ = mᵧ - nfree # obtain the system's symmetry vector for ω>0

    p = Q⁺ * (nᵀᵧ - mᵧ) # compute the vector p
    nᵀᵧꜛ⁰ = nᵀᵧ - Q * p # obtain the transverse symmetry vector for ω>0

    if nᵀᵧꜛ⁰ == mᵧꜛ⁰
        return all(pᵢ -> pᵢ ≈ round(pᵢ), p), p

    else
        return false, []
    end
end

function find_all_band_representations(
    vᵀ::AbstractSymmetryVector,
    long_modes::Vector{Vector{Int64}},
    d::Vector{Int64},
    brs::Collection{<:NewBandRep})

    # erase Γ from the high-symmetry points
    brs´ = prune_klab_irreps(brs, "Γ")
    vᵀ´ = prune_klab_irreps(vᵀ, "Γ")
    idxs = 1:length(first(brs´))


    # pick up only Γ from the high-symmetry points
    brsᵧ = pick_klab_irreps(brs, "Γ")
    vᵀᵧ = pick_klab_irreps(vᵀ, "Γ")

    # compute the fixed part and the free part of the physical ω=0 irreps at Γ
    lgirs = only(irreps(vᵀᵧ))
    klabel(first(lgirs)) == "Γ" || error("input symmetry vector to `physical` may only 
                                            reference Γ-contents")

    nfree, Q = physical_zero_frequency_gamma_irreps(
        lgirs;
        supergroup_constraints=true,
        force_fixed=true,
        lattice_reduce=true)

    # construct te vectors that will store the solutions
    phys_vec = Vector{Bool}[]
    p_vec = Vector{Vector{Float64}}[]
    solutions = Vector{Vector{Int64}}[]
    long_solutions = Vector{Int64}[]

    for i in eachindex(long_modes)
        nᴸ = long_modes[i]
        vᴸ´ = sum(brs´[nᴸ])
        vᵀ⁺ᴸ´ = vᵀ´ + vᴸ´
        μᵀ⁺ᴸ = occupation(vᵀ⁺ᴸ´)

        nᵀ⁺ᴸ = find_all_admissible_expansions(brs´, d, μᵀ⁺ᴸ, Vector(vᵀ⁺ᴸ´), idxs)

        if !isempty(nᵀ⁺ᴸ)
            check = [physical(vᵀᵧ, sum(brsᵧ[j]), sum(brsᵧ[nᴸ]), nfree, Q) for j in nᵀ⁺ᴸ]
            push!(solutions, nᵀ⁺ᴸ)
            push!(long_solutions, nᴸ)
            push!(phys_vec, [check[j][1] for j in eachindex(nᵀ⁺ᴸ)])
            push!(p_vec, [check[j][2] for j in eachindex(nᵀ⁺ᴸ)])
        end
    end
    return TightBindingCandidates(solutions, long_solutions, phys_vec, p_vec, brs)
end

function find_physical_band_representations(
    vᵀ::AbstractSymmetryVector,
    long_modes::Vector{Vector{Int64}},
    d::Vector{Int64}, brs::Collection{<:NewBandRep})
    all_solutions = find_all_band_representations(vᵀ, long_modes, d, brs)

    p_vec = Vector{Vector{Float64}}[]
    solutions = Vector{Vector{Int64}}[]
    long_solutions = Vector{Int64}[]

    for i in eachindex(all_solutions.phys)
        for j in eachindex(all_solutions.phys[i])
            if all_solutions.phys[i][j]
                push!(solutions, all_solutions.solutions[i])
                push!(long_solutions, all_solutions.long_modes[i])
                push!(p_vec, all_solutions.p[i])
            end
        end
    end

    return PhysicalTightBindingCandidates(solutions, long_solutions, p_vec, brs)
end