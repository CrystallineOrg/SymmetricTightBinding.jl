using Test
using SymmetricTightBinding
using SymmetricTightBinding: reciprocal_translation_phase, Wannier90TightBindingModel,
                             read_wannier90, fit_wannier90, _orbit_term_ranges,
                             _fit_hopping_orbit
using Crystalline
using LinearAlgebra
using StaticArrays: SVector

# deterministic, bounded test coefficients & momenta
_coefs(n; a = 0.73) = [0.3 * cospi(a * i) for i in 1:n]
_kpoints(D) = [[0.37 * cospi(0.71 * t * d) for d in 1:D] for t in 1:5]

# lattice-phase-only Bloch Hamiltonian (Convention 2): ∑_R e^{2πik·R} H(R)
function _conv2(w90tbm, k)
    return sum(zip(w90tbm.translations, w90tbm.hopping_matrices)) do (R, H)
        cispi(2 * dot(k, R)) * H
    end
end

# minimal `_hr.dat` reader: `(translations, hopping_matrices)` in order of appearance
function _read_hr(path)
    lines = readlines(path)
    N = parse(Int, lines[2])
    Nᴿ = parse(Int, lines[3])
    l = 4 + cld(Nᴿ, 15) # skip header and degeneracy weights
    translations = SVector{3, Int}[]
    hopping_matrices = Matrix{ComplexF64}[]
    for _ in 1:Nᴿ
        push!(translations, SVector{3, Int}(parse.(Int, split(lines[l])[1:3])))
        H = zeros(ComplexF64, N, N)
        for _ in 1:N^2
            xs = split(lines[l])
            i, j = parse.(Int, xs[4:5])
            H[i, j] = complex(parse(Float64, xs[6]), parse(Float64, xs[7]))
            l += 1
        end
        push!(hopping_matrices, H)
    end
    return translations, hopping_matrices
end

# check the `_hr.dat` layout of `tbm(cs)` against `tbm(cs)` itself, for two sets of `cs`
function test_w90_consistency(tbm)
    for a in (0.73, 0.41)
        ptbm = tbm(_coefs(length(tbm); a))
        w90tbm = Wannier90TightBindingModel(ptbm)
        @test orbital_positions(w90tbm) == orbital_positions(ptbm)
        @test orbital_count(w90tbm) == orbital_count(ptbm)
        @test allunique(w90tbm.translations)
        @test issorted(w90tbm.translations)
        for k in _kpoints(dim(ptbm))
            H = Matrix(ptbm(k))
            @test w90tbm(k) ≈ H
            Θ = reciprocal_translation_phase(orbital_positions(ptbm), k)
            @test _conv2(w90tbm, k) ≈ Θ' * H * Θ
        end
        if hermiticity(ptbm) == HERMITIAN # H(-R) = H(R)†
            Hᴿs = Dict(zip(w90tbm.translations, w90tbm.hopping_matrices))
            @test all(Hᴿs[-R] ≈ H' for (R, H) in Hᴿs)
        end
    end
end

@testset "Wannier90TightBindingModel" begin
    @testset "graphene (plane group 17)" begin
        brs = bandreps(17, Val(2))
        tbm = tb_hamiltonian(@composite(brs[5]), [[0, 0], [1, 0]])
        test_w90_consistency(tbm)
    end

    @testset "3D, incl. centered space group" begin
        brs = bandreps(16, Val(3))
        tbm = tb_hamiltonian(@composite(brs[1] + brs[3]), [[0, 0, 0], [1, 0, 0]])
        test_w90_consistency(tbm)

        brs = bandreps(227, Val(3)) # face-centered; positions in the primitive basis
        tbm = tb_hamiltonian(@composite(brs[1] + brs[2]), [[0, 0, 0]])
        test_w90_consistency(tbm)
    end

    @testset "1D: non-Hermitian & composite" begin
        brs = bandreps(2, Val(1))
        cbr = @composite brs[1] + brs[3]
        Rs = [[0], [1]]
        tbm_nh = tb_hamiltonian(cbr, Rs, Val(NONHERMITIAN))
        test_w90_consistency(tbm_nh)
        tbm_a = tb_hamiltonian(cbr, Rs, Val(ANTIHERMITIAN))
        test_w90_consistency(tbm_a)
        ctbm = tb_hamiltonian(cbr, Rs, Val(HERMITIAN)) + tbm_a
        test_w90_consistency(ctbm)
    end

    @testset "spinful" begin
        brs = bandreps(16; spinful = true, timereversal = false)
        tbm = tb_hamiltonian(@composite(brs[1] + brs[2]), [[0, 0, 0]])
        test_w90_consistency(tbm)
    end
end

@testset "write_wannier90" begin
    mktempdir() do dir
        # 2D: graphene, with lattice ⇒ also `_centres.xyz` & `.win`
        brs = bandreps(17, Val(2))
        tbm = tb_hamiltonian(@composite(brs[5]), [[0, 0], [1, 0]])
        ptbm = tbm(_coefs(length(tbm)))
        lattice = DirectBasis([1.0, 0.0], [-0.5, √3 / 2])
        prefix = joinpath(dir, "graphene")
        paths = write_wannier90(prefix, ptbm; lattice)
        @test paths == (; path_hr = prefix * "_hr.dat",
                          path_centres = prefix * "_centres.xyz",
                          path_win = prefix * ".win")

        # `_hr.dat` round-trip (R padded to 3D; elements at 12 decimals)
        w90tbm = Wannier90TightBindingModel(ptbm)
        translations, hopping_matrices = _read_hr(prefix * "_hr.dat")
        @test translations == [SVector(R[1], R[2], 0) for R in w90tbm.translations]
        @test all(isapprox.(hopping_matrices, w90tbm.hopping_matrices; atol = 1e-11))
        @test all(==("1"), split(readlines(prefix * "_hr.dat")[4])) # degeneracy weights

        # header: timestamped as by Wannier90 (e.g., "… on 04Nov2026 at 10:12:33"), with the
        # same header in all three files
        header = readlines(prefix * "_hr.dat")[1]
        timestamp = r"\d{2}[A-Z][a-z]{2}\d{4} at \d{2}:\d{2}:\d{2}$"
        @test occursin(r"^Created by SymmetricTightBinding\.jl on " * timestamp, header)
        @test readlines(prefix * "_centres.xyz")[2] == "Wannier centres. " * header
        @test readlines(prefix * ".win")[1] == "! " * header

        # centres: Cartesian positions, embedded in the xy-plane
        xyz = readlines(prefix * "_centres.xyz")
        @test parse(Int, xyz[1]) == 2
        for (line, r) in zip(xyz[3:end], orbital_positions(ptbm))
            tag, x, y, z = split(line)
            @test tag == "X"
            @test parse.(Float64, [x, y, z]) ≈ [cartesianize(r, lattice)..., 0.0] atol=1e-9
        end

        # `.win`: lattice vectors as rows, padded by a unit vector along z
        win = readlines(prefix * ".win")
        i = findfirst(==("begin unit_cell_cart"), win)
        @test win[i+1] == "ang"
        lattice_rows = [parse.(Float64, split(win[i+1+n])) for n in 1:3]
        @test lattice_rows ≈ [[1.0, 0.0, 0.0], [-0.5, √3 / 2, 0.0], [0.0, 0.0, 1.0]]
        @test win[i+5] == "end unit_cell_cart"

        # without a lattice, only `_hr.dat` is written; input may be a
        # `Wannier90TightBindingModel`
        prefix′ = joinpath(dir, "graphene_rs")
        @test write_wannier90(prefix′, w90tbm) ==
              (; path_hr = prefix′ * "_hr.dat", path_centres = "", path_win = "")
        @test readlines(prefix′ * "_hr.dat")[2:end] == readlines(prefix * "_hr.dat")[2:end]

        # dimension mismatch between model and lattice
        lattice³ = DirectBasis([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0])
        @test_throws TypeError write_wannier90(prefix, ptbm; lattice = lattice³)
    end
end

@testset "read_wannier90 & fit_wannier90" begin
    dir = mktempdir()

    # round-trip: `write_wannier90` → `read_wannier90` → `fit_wannier90` recovers `cs`
    brs², brs³ = bandreps(17, Val(2)), bandreps(227, Val(3))
    for (cbr, Rs) in ((@composite(brs²[5]), [[0, 0], [1, 0], [1, 1]]),  # graphene
                      (@composite(brs³[7]), [[0, 0, 0], [1, 0, 0]]))    # (16c|A₁g)
        tbm = tb_hamiltonian(cbr, Rs)
        ptbm = tbm(_coefs(length(tbm)))
        prefix = joinpath(dir, "roundtrip")
        write_wannier90(prefix, ptbm)
        w90tbm = read_wannier90(prefix, orbital_positions(ptbm))
        for k in _kpoints(dim(ptbm))
            @test w90tbm(k) ≈ ptbm(k) atol = 1e-10
        end
        @test fit_wannier90(tbm, w90tbm).cs ≈ ptbm.cs atol = 1e-10

        # a symmetry-breaking perturbation is projected out: the fit is an orthogonal
        # projection, so fitting the fitted model changes nothing
        perturbed = Wannier90TightBindingModel(
            w90tbm.translations,
            [H + 0.01 * cospi.(0.3 * (1:size(H, 1)) .+ 0.7 * (1:size(H, 2))') for H in w90tbm.hopping_matrices],
            w90tbm.positions)
        ptbm′ = fit_wannier90(tbm, perturbed)
        @test fit_wannier90(tbm, Wannier90TightBindingModel(ptbm′)).cs ≈ ptbm′.cs

        # the per-orbit problems partition the terms, & the round-trip fits are exact
        orbit_ranges = _orbit_term_ranges(tbm)
        @test reduce(vcat, collect.(orbit_ranges)) == 1:length(tbm)
        w90_hoppings = Dict(zip(w90tbm.translations, w90tbm.hopping_matrices))
        @test all(orbit_ranges) do terms
            orbit_fit = _fit_hopping_orbit(tbm, terms, w90_hoppings)
            orbit_fit.residual < 1e-10 * max(1, orbit_fit.w90_norm)
        end
    end

    # anti- & non-Hermitian models (incl. hermiticity-related blocks of the former)
    brs¹ = bandreps(2, Val(1))
    cbr¹ = @composite brs¹[1] + brs¹[3]
    for S in (ANTIHERMITIAN, NONHERMITIAN)
        tbm = tb_hamiltonian(cbr¹, [[0], [1]], S)
        ptbm = tbm(_coefs(length(tbm)))
        @test fit_wannier90(tbm, Wannier90TightBindingModel(ptbm)).cs ≈ ptbm.cs
    end

    # degeneracy weights & `_wsvec.dat` minimal-image shifts (1D model, embedded in 3D)
    prefix = joinpath(dir, "wsvec")
    write(prefix * "_hr.dat", """
        header
        1
        2
        1 2
        0 0 0 1 1  1.0 0.0
        1 0 0 1 1  0.5 0.0
        """)
    positions = [DirectPoint{1}([0.0])]
    w90tbm = read_wannier90(prefix, positions) # without `_wsvec.dat`
    @test w90tbm.translations == [SVector(0), SVector(1)]
    @test w90tbm.hopping_matrices ≈ [fill(1.0, 1, 1), fill(0.25, 1, 1)]
    write(prefix * "_wsvec.dat", """
        ## header
        0 0 0 1 1
        1
        0 0 0
        1 0 0 1 1
        2
        0 0 0
        -2 0 0
        """)
    w90tbm = read_wannier90(prefix, positions)
    @test w90tbm.translations == [SVector(-1), SVector(0), SVector(1)]
    @test w90tbm.hopping_matrices ≈ [fill(0.125, 1, 1), fill(1.0, 1, 1), fill(0.125, 1, 1)]
    @test_throws DimensionMismatch read_wannier90(prefix, [positions; positions])

    # malformed files: a missing `_wsvec.dat` entry; a truncated `_hr.dat`; a block of
    # records with a duplicate `(i, j)`; a block of records mixing translations
    write(prefix * "_wsvec.dat", "## header\n0 0 0 1 1\n1\n0 0 0\n")
    @test_throws ErrorException read_wannier90(prefix, positions)
    rm(prefix * "_wsvec.dat")
    write(prefix * "_hr.dat", "header\n1\n2\n1 2\n0 0 0 1 1  1.0 0.0\n")
    @test_throws ErrorException read_wannier90(prefix, positions)
    positions² = [DirectPoint{1}([0.0]), DirectPoint{1}([0.5])]
    write(prefix * "_hr.dat", """
        header
        2
        1
        1
        0 0 0 1 1  1.0 0.0
        0 0 0 2 1  0.5 0.0
        0 0 0 1 1  0.5 0.0
        0 0 0 2 2  1.0 0.0
        """)
    @test_throws ErrorException read_wannier90(prefix, positions²)
    write(prefix * "_hr.dat", """
        header
        2
        2
        1 1
        0 0 0 1 1  1.0 0.0
        0 0 0 2 1  0.5 0.0
        1 0 0 1 2  0.5 0.0
        0 0 0 2 2  1.0 0.0
        """)
    @test_throws ErrorException read_wannier90(prefix, positions²)
end
