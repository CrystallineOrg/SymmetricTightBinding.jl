using Pkg
Pkg.activate(@__DIR__)

# Design, rather than fitting: find hopping amplitudes that make the Γ₄⁺ triplet of a space
# group 224 (Pn-3m) model an *energy-isolated* three-fold degeneracy, i.e., one whose energy
# no other band attains anywhere in the Brillouin zone. Cf. `docs/src/design.md`.

using Crystalline, SymmetricTightBinding, Optim, Random, Brillouin, GLMakie

sgnum, D = 224, 3
brs = calc_bandreps(sgnum, Val(D))
cbr = @composite brs[18] + brs[23]           # (4b|A₂g) ⊕ (2a|A₂): 6 bands
tbm = tb_hamiltonian(cbr, [[0, 0, 0], [1, 0, 0]])

# the Γ-point content is 2Γ₂⁺ + Γ₄⁺ + Γ₁⁻: a single 3-dimensional irrep and three singlets,
# so isolating Γ₄⁺ means pushing all three singlets to one side of the triplet's energy
println(SymmetryVector(cbr))

# NB: `Nk` divisible by 6, so that the mesh samples the ½- *and* ⅓-type high-symmetry
# points; and many starts — design searches are much harder than fits (50 finds nothing)
Random.seed!(1)
ptbm = isolate_irrep(tbm, "Γ₄⁺"; Nk = 6, max_multistarts = 200, verbose = true)

# verify on a mesh 4× finer than the one optimized over
report = isolation_report(ptbm, "Γ₄⁺"; Nk = 24)
println(report) # → 2 of the triplet's bands below Et, 1 above; i.e., a filling of 4

# refining against a denser mesh is much cheaper warm-started from the result above
Random.seed!(2)
ptbm_fine = isolate_irrep(tbm, "Γ₄⁺"; Nk = 12, max_multistarts = 20, init = ptbm)
println(isolation_report(ptbm_fine, "Γ₄⁺"; Nk = 24))

# a `lasso` weight asks for the sparsest model that still isolates the multiplet (warm-
# starting from the result above, since switching the penalty on changes the landscape only
# gradually); it is a blunt instrument — sweep the weight rather than trusting one value
Random.seed!(1)
ptbm_sparse = isolate_irrep(tbm, "Γ₄⁺"; Nk = 6, max_multistarts = 30, lasso = 0.4,
                            init = ptbm)
println(isolation_report(ptbm_sparse, "Γ₄⁺"; Nk = 24))
println("nonzero amplitudes: ", count(>(1e-3), abs.(ptbm_sparse.cs)), " of ", length(tbm))

kpi = interpolate(irrfbz_path(sgnum, directbasis(sgnum, Val(D))), 100)
faxp = plot(
    kpi, spectrum(ptbm, kpi), fill(report.Et, length(kpi), 1);
    color = [:black, :crimson], linewidth = [3, 1], linestyle = [nothing, :dash],
    label = ["Bands", rich("E", subscript("t"))]
)
axislegend(; framevisible = false, position = :rb)
display(faxp)
