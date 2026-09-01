using Crystalline, SymmetricTightBinding, GLMakie, Brillouin, Optim, Random
using LinearAlgebra

sgnum = 224
brs = calc_bandreps(sgnum)
cbr = @composite brs[18] + brs[23]
tbm = tb_hamiltonian(cbr, [[0, 0, 0], [1, 0, 0]]) # on-site + nearest neighbour
SymmetryVector(cbr)

kp = irrfbz_path(224, directbasis(224))
kpi = interpolate(kp, 500);
kpi_coarse = interpolate(kp, 100);
#ks = filter(k->norm(k) > .1, kpi_coarse)


#Random.seed!(1)
#ptbm = isolate_irrep(tbm, "Γ₄⁺"; max_multistarts=10, verbose=true, ks=kpi_coarse)

# fit a reduced model, with just a subset of terms (found by ad hoc elimination/pruning)
tbm_r = tbm[[1, 4, 7, 12, 10]] # 10 is the least relevant one, but something is needed to break Gamma-R equifrequency
#tbm_r = tbm[[4, 12, 10]] # also gives a quite good 3-fold degeneracy at Gamma - and, intriguingly, also a good one 3-fold one at R
ptbm_r = isolate_irrep(tbm_r, "Γ₄⁺"; max_multistarts=10, verbose=true, ks=kpi_coarse)
plot(kpi, spectrum(ptbm_r, kpi); annotations=collect_irrep_annotations(ptbm_r))

# add some term that break symmetry
comp_tbm = subduced_complement(tbm, 111)

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

#=
new_c = -0.05
for n in 1:length(comp_tbm)
    tbm_broken = vcat(tbm_r, comp_tbm[n:n])
    ptbm_broken = ParameterizedTightBindingModel{3, HERMITIAN}(tbm_broken, vcat(ptbm_r.cs, new_c))

    faxp = plot(kpi, spectrum(ptbm_r, kpi), spectrum(ptbm_broken, kpi);
                color = [:gray80, "rgb(39,60,117)"])
    display(faxp)
    wait_for_key("model $n/$(length(comp_tbm))")
end
=#

#=
new_c = -0.02
for n in 1:length(comp_tbm)
    tbm_broken = vcat(tbm_r, comp_tbm[10:10], comp_tbm[n:n])
    ptbm_broken = ParameterizedTightBindingModel{3, HERMITIAN}(tbm_broken, vcat(ptbm_r.cs, 0.03, new_c))

    faxp = plot(kpi, spectrum(ptbm_r, kpi), spectrum(ptbm_broken, kpi);
                color = [:gray80, "rgb(39,60,117)"])
    display(faxp)
    wait_for_key("model $n/$(length(comp_tbm))")
end
=#

#=
new_c10 = 0.03
tbm_broken_simple = vcat(tbm_r, comp_tbm[[10]])
ptbm_broken_simple = ParameterizedTightBindingModel{3, HERMITIAN}(tbm_broken_simple, vcat(ptbm_r.cs, new_c10))
n = 12
for new_cn in range(0, -0.1, 11)
    tbm_broken = vcat(tbm_r, comp_tbm[[10, n]]);
    ptbm_broken = ParameterizedTightBindingModel{3, HERMITIAN}(tbm_broken, vcat(ptbm_r.cs, new_c10, new_cn));
    faxp = plot(kpi, spectrum(ptbm_r, kpi), spectrum(ptbm_broken_simple, kpi), spectrum(ptbm_broken, kpi);
                color = [:gray80, :red, "rgb(39,60,117)"], ylims = (-0.1, 1.2))
    display(faxp)
    wait_for_key("   new_c12 = $new_cn")
end
=#

# Final "symmetry-breaking" design with nice isolated Weyl points (only optimal for the 
# `tbm_r = tbm[[1, 4, 7, 12, 10]]` case above)

c10 = 0.03
c12 = -0.06

tbm_broken = vcat(tbm_r, comp_tbm[[10, 12]]);
ptbm_broken = ParameterizedTightBindingModel{3, HERMITIAN}(tbm_broken, vcat(ptbm_r.cs, c10, c12));
ptbm_broken_simple = ParameterizedTightBindingModel{3, HERMITIAN}(ptbm_broken.tbm[1:end-1], ptbm_broken.cs[1:end-1])
faxp = plot(
    kpi, 
    spectrum(ptbm_r, kpi), #=spectrum(ptbm_broken_simple, kpi), =#spectrum(ptbm_broken, kpi);
    color = [:gray80, #=:crimson, =# "rgb(39,60,117)"],
    ylims = (-0.1, 1.2)
)

Es = range(.5, .75, 200)
dos = densityofstates(ptbm_r, Es; Nk = 200)
dos_broken = densityofstates(ptbm_broken, Es; Nk = 200)
lines(Es, dos; color=:gray80)
lines!(Es, dos_broken; color="rgb(39,60,117)")
