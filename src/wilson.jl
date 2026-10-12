using LinearAlgebra: I, dot, eigvals!, svd!, mul!

#=
using Crystalline, SymmetricTightBinding, GLMakie
sgnum = 13
brs = bandreps(sgnum, Val(2); timereversal = false)
brs[4], brs[1] = brs[1], brs[4];
cbr = @composite brs[1] + brs[4] # (1b|A) + (1c|A)
tbm = tb_hamiltonian(cbr, [[0,0], [0,1]])
 haldane_model(t₁, m, t₂, ϕ) = tbm([m, t₂*cos(ϕ), -t₂*sin(ϕ), -m, t₂*cos(ϕ), t₂*sin(ϕ), t₁, 0])
ptbm = haldane_model(1.0, 0.0, 0.1, π/2) # inversion symmetric

k12s = range(-0.5, 0.5, 50)
ksv = ((ReciprocalPoint(k1, k2) for k2 in k12s) for k1 in k12s)

#cd("../../../../../mnt/c/Users/tchr/OneDrive - Danmarks Tekniske Universitet/Files/Projects/photonic-topo-sym/julia/SymmetricTightBinding")
#const STB = SymmetricTightBinding
#@eval STB include("src/wilson.jl")
#@eval STB export wilsonloop

ps = wilsonloop.(Ref(ptbm), ksv, Ref(1:2))
lines(kxys, ps)
=#
function wilsonloop(
    ptbm :: ParameterizedTightBindingModel{D},
    ks, # iterable of D-length `AbstractVector{<:Real}`s
    bands :: AbstractVector{<:Integer};
    scratch_spaces :: NTuple{3, Matrix{ComplexF64}} = begin # *three* scratch matrices
        # re-uses the default buffer from `ptbm` for one `scratch_spaces[3]`: might as well
        wilsonloop_scratch_spaces(ptbm)
    end
) where {D}
    hermiticity(ptbm) === HERMITIAN || error("anti- & non-hermitian models not yet supported")
    N = orbital_count(ptbm)
    scratch1, scratch2, us₁_copy = scratch_spaces

    # validate input
    if !all(k -> eltype(k) <: Real && length(k) == D, ks)
        error(lazy"input `ks` must be iterable of real $D-vectors")
    end
    G = last(ks) - first(ks)
    if !all(Gᵢ -> isapprox(round(Gᵢ), Gᵢ; atol=VEC_CMP_ATOL), G)
        error(lazy"input `ks` does not form a loop (`ks[end] - ks[begin] = $G`)")
    end
    if !all(b-> 1≤b≤N, bands)
        error("input `bands` contains indices beyond the orbital dimension of `tbm`")
    end
    if !(size(scratch1) == size(scratch2) == size(us₁_copy) == (N, N))
        error("invalid input `scratch_spaces`: must be equal size and compatible with `ptbm`")
    end
    
    # setup for stepping through the `ks` loop
    bloch_phase = Val(false)
    k₁, ks_rest = Iterators.peel(ks) # `ks_rest = ks[2:end]`
    _, usᵢ₋₁ = solve(ptbm, k₁; scratch=scratch1, bloch_phase)
    copy!(us₁_copy, usᵢ₋₁)
    Nk = length(ks)
    
    # step through `ks` loop, building up the link-matrix `U` step by step
    U = Matrix{ComplexF64}(I(length(bands)))
    scratchᵢ = scratch2
    for kᵢ in Iterators.take(ks_rest, Nk-2) # ≡`ks[2:end-1]`, but work for iterators too
        _, usᵢ = solve(ptbm, kᵢ; scratch=scratchᵢ, bloch_phase) # get eigenvectors at kᵢ
        # next "step" for the link matrices
        U = update_link_matrices!(U, bands, usᵢ₋₁, usᵢ)

        # swap out scratch-space and eigenvector variables, cf. mutation in `solve`
        scratchᵢ = scratchᵢ === scratch1 ? scratch2 : scratch1
        usᵢ₋₁ = usᵢ
    end

    # last k-point equals the first, mod G
    Θᴳ = reciprocal_translation_phase(orbital_positions(ptbm), G)
    mul!(scratchᵢ, Θᴳ, us₁_copy); usᵢ = scratchᵢ # multiply by exp(-iG⋅r)
    U = update_link_matrices!(U, bands, usᵢ₋₁, usᵢ)

    # get phases from link matrices
    phases = sort!(angle.(eigvals!(unitarize!(U))))

    return phases
end

function wilsonloop(ptbm::ParameterizedTightBindingModel, ks, bands::Integer; kws...)
    return only(wilsonloop(ptbm, ks, [bands]; kws...))
end

function update_link_matrices!(U, bands, usᵢ₋₁, usᵢ)
    if length(bands) == 1 # small optimization for Abelian case
        band = first(bands)
        U[1] *= dot((@view usᵢ[:, band]), (@view usᵢ₋₁[:, band]))
    else # non-Abelian case
        U .= (@view usᵢ[:, bands])' * (@view usᵢ₋₁[:, bands]) * U
        # NB: a slightly more efficient thing could be done by passing two additional 
        #     `U`-like work-matrices here and writing the temporaries into them - but it's
        #     not worthwhile in the scheme of things, as ~80% of time is spent in `solve`
    end
    return U
end
function unitarize!(A)
    F = svd!(A) # uses `A` as buffer; write to `A` again below
    return mul!(A, F.U, F.Vt)
end

function wilsonloop_scratch_spaces(ptbm)
    scratch3 = ptbm.scratch # reuse the scratch space from `ptbm`; might as well
    scratch1, scratch2 = similar(scratch3), similar(scratch3)
    return (scratch1, scratch2, scratch3)
end