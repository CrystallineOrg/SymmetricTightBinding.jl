# Shared definitions for the analysis-environment scripts (`01_setup.jl`, `03_fit.jl`)

using Crystalline, SymmetricTightBinding, LinearAlgebra
using SymmetricTightBinding: site_induced_sgrep_excl_phase, primitivized_generators,
                             Wannier90TightBindingModel

const DOJO_LDA = "dojo.nc.sr.lda.v0_4_1.standard.upf" # PseudoDojo, via PseudoPotentialData.jl

# ---------------------------------------------------------------------------------------- #
# Cases
#
# - `atoms`: (species, Wyckoff position label) pairs
# - `brlabels`: the EBRs of the Wannierized bands; `angular[brlabel]`: the angular part
#   `(l, m)` of the initial projection, per partner function of the site irrep (see
#   `projection_angular_parts`)
# - `dft["reference"]`: energy zero; "vbm" (valence-band maximum; insulators) or "fermi"
# - `wannier`: either `bands` (an isolated set of bands, by band index; the outer & frozen
#   windows are then placed in the gaps around them, if needed) or `n_bands` & windows
#   relative to the energy zero (disentanglement); and `num_iter` (0: projected Wannier
#   functions)

const CASES = Dict(
    "si-valence" => (;
        description = "Si: the 4 valence bands, as bond-centred s-like Wannier functions",
        sgnum = 227, dim = 3, a = 5.431, # Å, experimental lattice constant
        atoms = [("Si", "8a")],
        brlabels = ["(16c|A₁g)"],
        angular = Dict("(16c|A₁g)" => [(0, 0)]),
        dft = Dict("Ecut" => 16, "kgrid" => [8, 8, 8], "pseudo_family" => DOJO_LDA,
                   "n_bands" => 8, "kline_density" => 40, "reference" => "vbm"),
        wannier = Dict("bands" => [1, 4], "num_iter" => 1000),
    ),
    "hbn-pz" => (;
        description = "hBN monolayer: the two π bands, as pz Wannier functions on N and B",
        sgnum = 14, dim = 2, a = 2.504, vacuum = 12.0, # Å; plane group p3m1
        atoms = [("N", "1c"), ("B", "1b")],
        brlabels = ["(1c|A₁)", "(1b|A₁)"], # pz: invariant under the in-plane operations
        angular = Dict("(1c|A₁)" => [(1, 0)], "(1b|A₁)" => [(1, 0)]),
        dft = Dict("Ecut" => 30, "kgrid" => [12, 12, 1], "pseudo_family" => DOJO_LDA,
                   "n_bands" => 16, "kline_density" => 40, "reference" => "vbm"),
        # the π bands are entangled in energy with the σ bands, but decoupled from them by
        # the mirror z → -z (pz is odd, σ is even), which the projections respect; frozen
        # window: π & π* near K and M only (cf. `dft_bands.log`)
        wannier = Dict("n_bands" => 16, "num_iter" => 1000, "dis_num_iter" => 1000,
                       "dis_froz_min" => -1.0, "dis_froz_max" => 4.7),
    ),
    "graphene-pz" => (;
        description = "graphene: the two π bands, as pz Wannier functions on C",
        sgnum = 17, dim = 2, a = 2.461, vacuum = 12.0, # Å; plane group p6mm
        atoms = [("C", "2b")],
        brlabels = ["(2b|A₁)"], # pz: invariant under the in-plane operations
        angular = Dict("(2b|A₁)" => [(1, 0)]),
        dft = Dict("Ecut" => 30, "kgrid" => [12, 12, 1], "pseudo_family" => DOJO_LDA,
                   "temperature" => 0.002, # Ha; Gaussian smearing (semimetal)
                   "n_bands" => 16, "kline_density" => 40, "reference" => "fermi"),
        # as for hBN; frozen window: from just above the top of the σ bands (E_F − 3.02 eV,
        # at Γ) to just below the lowest other unoccupied state (E_F + 3.35 eV, at Γ), i.e.,
        # π & π* around K and M (cf. `dft_bands.log`)
        wannier = Dict("n_bands" => 16, "num_iter" => 1000, "dis_num_iter" => 1000,
                       "dis_froz_min" => -2.9, "dis_froz_max" => 3.0),
    ),
    "srvo3-t2g" => (;
        description = "SrVO₃: the 3 t₂g bands at the Fermi level, as V t₂g Wannier functions",
        sgnum = 221, dim = 3, a = 3.842, # Å; cubic perovskite
        atoms = [("Sr", "1b"), ("V", "1a"), ("O", "3d")],
        brlabels = ["(1a|T₂g)"],
        angular = Dict("(1a|T₂g)" => [(2, -1), (2, 1), (2, -2)]), # d_yz, d_xz, d_xy
        dft = Dict("Ecut" => 30, "kgrid" => [8, 8, 8], "pseudo_family" => DOJO_LDA,
                   "temperature" => 0.002, # Ha; Gaussian smearing (metal)
                   "n_bands" => 26, "kline_density" => 40, "reference" => "fermi"),
        # the t₂g bands (21–23) overlap in energy with the e_g bands above (but not with
        # the O 2p bands below; cf. `dft.log`): disentanglement, with a frozen window below
        # the bottom of the e_g bands (at E_F + 1.19 eV)
        wannier = Dict("n_bands" => 26, "num_iter" => 1000, "dis_num_iter" => 1000,
                       "dis_win_min" => -1.6, "dis_win_max" => 3.0,
                       "dis_froz_min" => -1.6, "dis_froz_max" => 1.15),
    ),
)

# ---------------------------------------------------------------------------------------- #
# Band representations & sites

# find a band representation by its label, e.g., "(16c|A₁g)"
function findbr(brs, label::AbstractString)
    idx = findall(br -> string(br) == label, brs)
    length(idx) == 1 || error("no unique BR with label $label")
    return only(idx)
end

# composite BR from a list of EBR labels (each included once)
function composite_from_labels(brs, labels)
    coefs = zeros(Int, length(brs))
    for l in labels
        coefs[findbr(brs, l)] += 1
    end
    return CompositeBandRep(coefs, brs)
end

# the sites of the Wyckoff position labelled `wp_label` (e.g., "8a"), in the primitive
# basis (as the orbital positions of any BR at that Wyckoff position)
function wyckoff_sites(brs, wp_label::AbstractString)
    idx = findfirst(br -> Crystalline.label(position(br)) == wp_label, brs)
    isnothing(idx) && error("no BR at Wyckoff position $wp_label")
    return unique(orbital_positions(brs[idx]))
end

# initial Wannier90 projections, as `(l, m)` per orbital, in SymmetricTightBinding's
# orbital order: BRs in the order of `cbr.brs` and, within each BR, partner-function-major
# (all sites of partner function 1, then all sites of partner function 2, ...).
# `angular[label]` lists `(l, m)` for each partner function of the site irrep of `label`,
# with DFTK's real harmonics: `(0, 0)` = s; `(1, 1)`, `(1, -1)`, `(1, 0)` = pₓ, p_y, p_z;
# `(2, -2)`, `(2, -1)`, `(2, 1)` = d_xy, d_yz, d_xz; `(2, 0)`, `(2, 2)` = d_z², d_x²-y²
function projection_angular_parts(cbr, angular::Dict)
    lms = Tuple{Int, Int}[]
    for (br, c) in zip(cbr.brs, cbr.coefs)
        iszero(c) && continue
        lms_br = angular[string(br)]
        nsites = occupation(br) ÷ length(lms_br)
        for _ in 1:c, lm in lms_br, _ in 1:nsites
            push!(lms, lm)
        end
    end
    return lms
end

# ---------------------------------------------------------------------------------------- #
# Aligning the local orbital basis of the Wannier functions with that of the tight-binding
# model: the Wannier functions are (approximately) the real-harmonic orbitals of the
# projections, `(l, m)` per orbital, while SymmetricTightBinding's orbitals transform under
# the site-symmetry group via Crystalline's site irreps (and are related across the sites
# of a Wyckoff position via coset representatives). The two bases are related by a
# site-local unitary `U` (|stbₐ⟩ = ∑_b |w90_b⟩ U[b,a]), which intertwines the two
# representations of the space group, i.e., `ρ_w90(g) U = U ρ_stb(g)`; by Schur's lemma, `U`
# is unique up to a phase per site irrep.

# DFTK's real solid harmonics (`DFTK.solid_harmonic_real`), for l ≤ 2
function solid_harmonic_real(l, m, r)
    x, y, z = r
    l == 0 && return sqrt(1 / 4π)
    if l == 1
        return sqrt(3 / 4π) * (m == -1 ? y : m == 0 ? z : x)
    elseif l == 2
        m == -2 && return sqrt(15 / 4π) * x * y
        m == -1 && return sqrt(15 / 4π) * y * z
        m ==  0 && return sqrt(5 / 16π) * (2z^2 - x^2 - y^2)
        m ==  1 && return sqrt(15 / 4π) * x * z
        m ==  2 && return sqrt(15 / 16π) * (x^2 - y^2)
    end
    error("l > 2 not implemented")
end

# representation `Dˡ` of a Cartesian rotation `W` (3×3) on the real harmonics of degree `l`:
# (g Yₘ)(r) = Yₘ(W⁻¹r) = ∑ₘ′ Yₘ′(r) Dˡ[m′,m], found by sampling (exact: the l-shell is
# closed under rotations); rows/columns are ordered as m = -l, …, l
function realharmonic_rotation(l, W)
    rs = [[cospi(0.37i) * sinpi(0.21i + 0.1), sinpi(0.37i) * sinpi(0.21i + 0.1),
           cospi(0.21i + 0.1)] for i in 1:4(2l+1)] # generic sample points
    Y  = [solid_harmonic_real(l, m, r) for r in rs, m in -l:l]
    Y′ = [solid_harmonic_real(l, m, W' * r) for r in rs, m in -l:l] # W⁻¹ = Wᵀ
    return Y \ Y′
end

# Cartesian, 3D form of a (primitive-basis) point operation `W` on a `D`-dimensional lattice
# (columns of `lattice`); for `D < 3`, the extra dimensions are left invariant
function cartesian_rotation(W, lattice)
    D = size(lattice, 1)
    W³ = Matrix{Float64}(I, 3, 3)
    W³[1:D, 1:D] .= lattice * W / lattice
    return W³
end

# representation of a (primitive-basis) space group operation `g` on the real-harmonic
# orbitals `lms` at `positions` (excluding lattice-translation phases)
function realharmonic_rep(g, positions, lms, lattice)
    W = rotation(g)
    W³ = cartesian_rotation(W, lattice)
    Dˡs = Dict(l => realharmonic_rotation(l, W³) for l in unique(first.(lms)))
    N = length(positions)
    ρ = zeros(N, N)
    for a in 1:N
        ga = W * positions[a] + Crystalline.translation(g) # site that `g` maps site `a` to
        for b in 1:N
            d = ga - positions[b]
            (lms[a][1] == lms[b][1] && isapprox(d, round.(d); atol = 1e-8)) || continue
            l = lms[a][1]
            ρ[b, a] = Dˡs[l][lms[b][2]+l+1, lms[a][2]+l+1]
        end
    end
    return ρ
end

function basis_alignment(cbr, lms, lattice)
    positions = orbital_positions(cbr)
    N = length(positions)
    gens = primitivized_generators(cbr.brs[findfirst(!iszero, cbr.coefs)])
    # solve `ρ_w90(g) U - U ρ_stb(g) = 0` for site-local `U` (`U[b,a] = 0` if rₐ ≠ r_b)
    site_local = vec([isapprox(positions[b], positions[a]) for b in 1:N, a in 1:N])
    M = reduce(vcat, [kron(I(N), realharmonic_rep(g, positions, lms, lattice)) -
                      kron(transpose(site_induced_sgrep_excl_phase(cbr, g)), I(N))
                      for g in gens])[:, site_local]
    F = svd(M)
    nullity = count(<(1e-10), F.S)
    nullity > 0 || error("no intertwiner: Wannier90 projections and BR are incompatible")
    X = zeros(ComplexF64, N, N) # a generic intertwiner: combination of the null vectors
    for i in 1:nullity
        X[site_local] .+= F.V[:, end-i+1] ./ (i + 0.618)
    end
    U = X / sqrt(Hermitian(X' * X)) # unitary (polar) factor, still an intertwiner
    # fix the phase freedom per BR (i.e., per site irrep): make the first nonzero element of
    # the BR's first column real & positive
    offset = 0
    for (br, c) in zip(cbr.brs, cbr.coefs), _ in 1:c
        cols = offset+1:offset+occupation(br)
        u₁ = U[findfirst(x -> abs(x) > 1e-8, U[:, first(cols)]), first(cols)]
        U[:, cols] ./= u₁ / abs(u₁)
        offset += occupation(br)
    end
    maximum(abs, imag(U)) < 1e-10 && (U = real(U))
    return U
end

# hopping matrices of `w90tbm` in the basis of the tight-binding orbitals: `U' H(R) U`
function align_basis(w90tbm, U)
    return Wannier90TightBindingModel(w90tbm.translations,
                                      [Matrix{ComplexF64}(U' * H * U) for H in w90tbm.hopping_matrices],
                                      w90tbm.positions)
end
