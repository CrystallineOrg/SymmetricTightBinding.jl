using LinearAlgebra: I, dot, eigvals!, svd!, mul!

"""
    wilsonloop(
        ptbm :: ParameterizedTightBindingModel,
        ks,
        bands :: Union{Integer, AbstractVector{<:Integer}}
    ) --> Float64 or Vector{Float64}

Return the Wilson loop spectrum, or Berry phases, in (−π, π] of `ptbm` along a discretized
**k**-space loop `ks` over the specified `bands` (a single band index or an abstract vector
of indices). The `bands` must be gapped from all other bands everywhere along `ks`.

The loop `ks` must be closed, i.e., `last(ks) - first(ks)` must be a (primitive) reciprocal
lattice vector (for contractible loops, zero), and `ks` must be an iterable of
`dim(ptbm)`-length real-valued abstract vectors or tuples.

## Output

- If `bands` is a single integer, the Wilson loop is Abelian and a single phase is returned.
- If `bands` is a vector of integers, the Wilson loop is non-Abelian and a sorted vector of
  phases is returned.

## Examples

**Abelian Berry phases:**
A Haldane-like model with a Chern number of -1:

```julia
sgnum = 13 # p3
brs = bandreps(sgnum, Val(2); timereversal = false)
cbr = @composite brs[1] + brs[4] # (1c|A) + (1b|A)

tbm = tb_hamiltonian(cbr, [[0,0], [1,0]])
ptbm = tbm([0.0, 0.1, -0.1, 0.0, 0.1, 0.1, 1.0, 0.0]) # Haldane-like model

k12s = range(-0.5, 0.5, 101)
ksv = (([k1, k2] for k2 in k12s) for k1 in k12s)
phases = [wilsonloop(ptbm, ks, 1) for ks in ksv]
```

If plotted with e.g., `GLMakie.plot(k12s, phases)` this shows a down-directed 2π-winding
of the associated (Abelian) Berry phase, corresponding to a Chern number of -1.

**Non-Abelian Berry phases:**
The Kane-Mele model is a classical example of a ℤ₂ topological insulator whose topology can
be diagnosed by the non-Abelian Berry phases of its degenerate valence band doublet.

```julia
sgnum = 17 # p6mm
brs = bandreps(sgnum, Val(2); timereversal = true, spinful = Val(true))
cbr = @composite brs[3] # (2b|E₁ˢ)

tbm = tb_hamiltonian(cbr, [[0,0], [1,0]]) # nearest- and next-nearest-neighbor terms
tbm = tbm[[3,6]] # keep only "graphene" nearest-neighbor and Kane-Mele spin-orbit terms
ptbm = tbm([1.0, 0.05]) # Kane-Mele model with 0.05 spin-orbit coupling

phases = [wilsonloop(ptbm, ks, 1:2) for ks in ksv] # Wilson loop spectra for valence doublet
```

If plotted (e.g., `GLMakie.series(k12s, stack(phases))`) this shows the signature relative
winding of the non-Abelian Berry phases of the Kane-Mele model, indicating a non-trivial ℤ₂
topological invariant.

# Extended help

### Implementation details

The Wilson loop operator ``\\mathbf{W}`` is computed as the product of link matrices along
the discretized **k**-space loop ``\\mathbf{k}_1, …, \\mathbf{k}_N``[^1][^2],

```math
\\mathbf{W} = \\mathbf{U}(N→N-1) \\mathbf{U}(N-1→N-2) ⋯ \\mathbf{U}(2→1),
```

with each link matrix ``\\mathbf{U}(i→i-1)`` evaluated as the overlap of the Bloch
eigenvectors between adjacent **k**-points:

```math
U_{bb'}(i→i-1) = ⟨u_b(\\mathbf{k}_i)|u_{b'}(\\mathbf{k}_{i-1})⟩,
```

with ``|u_b(\\mathbf{k})⟩`` denoting the Bloch eigenvector at ``\\mathbf{k}`` in band ``b``.
The loop is closed by fixing the gauge at ``\\mathbf{k}_N = \\mathbf{k}_1 + \\mathbf{G}`` to
``|u_b(\\mathbf{k}_N)⟩ = \\mathrm{e}^{-\\mathrm{i}\\mathbf{G}⋅\\mathbf{r}}
|u_b(\\mathbf{k}_1)⟩``, with ``\\mathbf{r}`` the orbital positions.

The associated Wilson loop spectrum (Berry phases) is computed as the polar angles of the
eigenvalues of ``\\mathbf{W}``, after projecting ``\\mathbf{W}`` onto the nearest unitary
matrix (via its singular value decomposition).

[^1]: Alexandradinata, Dai, & Bernevig, *Wilson-loop characterization of
    inversion-symmetric topological insulators*,
    [Phys. Rev. B **89**, 155114 (2014)](https://doi.org/10.1103/PhysRevB.89.155114).

[^2]: Gresch *et al.*, *Z2Pack: Numerical implementation of hybrid Wannier centers for
    identifying topological materials*,
    [Phys. Rev. B **95**, 075146 (2017)](https://doi.org/10.1103/PhysRevB.95.075146).
"""
function wilsonloop(
    ptbm :: ParameterizedTightBindingModel{D},
    ks, # iterable of D-length `AbstractVector{<:Real}`s
    bands :: AbstractVector{<:Integer};
) where {D}
    N = orbital_count(ptbm)

    # validate input
    hermiticity(ptbm) === HERMITIAN || error("anti- & non-hermitian models not yet supported")
    if !all(k -> eltype(k) <: Real && length(k) == D, ks)
        error(lazy"input `ks` must be iterable of real $D-vectors")
    end
    G = ReciprocalPoint{D}(last(ks)) - ReciprocalPoint{D}(first(ks)) # Reci… for `ks[i]::Tuple`
    if !all(Gᵢ -> isapprox(round(Gᵢ), Gᵢ; atol=VEC_CMP_ATOL), G)
        error(lazy"input `ks` does not form a loop (`ks[end] - ks[begin] = $G`)")
    end
    if !all(b-> 1≤b≤N, bands)
        error("input `bands` contains indices beyond the orbital dimension of `tbm`")
    end
    length(bands) == 0 && error("input `bands` must contain at least one band index")
    allunique(bands) || error("input `bands` must contain only unique band indices")
    length(ks) ≤ 1 && error("input `ks` must contain at least 2 k-points to form a loop")

    # 2 scratch matrices: none shared with `ptbm.scratch`, so threading is safe
    scratch1, scratch2 = similar(ptbm.scratch), similar(ptbm.scratch)

    # setup for stepping through the `ks` loop
    bloch_phase = Val(false)
    k₁, ks_rest = Iterators.peel(ks) # `ks_rest = ks[2:end]`
    _, usᵢ₋₁ = solve(ptbm, k₁; scratch=scratch1, bloch_phase)
    us₁_copy = copy(usᵢ₋₁) # needed for last k-point, after loop
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

    # last k-point equals the first, mod G; link via u_kᴺ(r) = exp(-iG⋅r)u_k¹(r)
    Θᴳ = reciprocal_translation_phase(orbital_positions(ptbm), G) # = exp(-iG⋅r)
    mul!(scratchᵢ, Θᴳ, us₁_copy); usᵢ = scratchᵢ # multiply by exp(-iG⋅r)
    U = update_link_matrices!(U, bands, usᵢ₋₁, usᵢ)

    # get phases from link matrices
    phases = sort!(angle.(eigvals!(unitarize!(U))))

    return phases
end

# return a scalar `Float64` phase if `bands` is a single integer
function wilsonloop(ptbm::ParameterizedTightBindingModel, ks, bands::Integer)
    return only(wilsonloop(ptbm, ks, bands:bands))
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

# the Vanderbilt trick; obtain a strictly unitary matrix from a near-unitary matrix
function unitarize!(A)
    F = svd!(A) # uses `A` as buffer; write to `A` again below
    return mul!(A, F.U, F.Vt)
end