"""
    spectrum(ptbm::ParameterizedTightBindingModel, ks)

Evaluate the spectrum, i.e., energies, of the tight-binding model `ptbm` over an iterable
of input momenta `ks`. 

The energies are returned as a matrix, with rows running over the bands and columns over
momenta.

## Example

As an example, we evaluating the band structure of graphene.
Below, we first construct and parameterize a tight-binding model for the the (2b|A₁) EBR in
plane group 17, corresponding to the highest-lying orbitals in graphene.
Next, we construct a path along high-symmetry directions of the Brillouin zone using
[Brillouin.jl](https://github.com/thchr/Brillouin.jl), calculate the spectrum across this
path; and finally, plot the band structure using Brillouin and PlotlyJS:

```julia-repl
julia> using Crystalline, TETB

julia> brs = calc_bandreps(17, Val(2)); coefs = zeros(Int, length(brs)); coefs[5] = 1; 

julia> cbr = CompositeBandRep(coefs, brs)
13-irrep CompositeBandRep{2}:
 (2b|A₁) (2 bands)

julia> ptbm = tb_hamiltonian(cbr, [zeros(Int, dim(cbr))])([0.0, 1.0]);

julia> using Brillouin, PlotlyJS

julia> kpi = interpolate(irrfbz_path(17, directbasis(17, Val(2))), 100);

julia> plot(kpi, spectrum(ptbm, kpi)')
```
"""
function spectrum(ptbm::ParameterizedTightBindingModel, ks)
    if !(eltype(ks) <: AbstractVector{<:Real})
        error("the elements of `ks` must subtype `AbstractVector{<:Real}`")
    end
    Es = Matrix{Float64}(undef, ptbm.tbm.N, length(ks))
    for (i, k) in enumerate(ks)
        es = spectrum(ptbm, k)
        @inbounds Es[:, i] .= es
    end
    return Es
end

"""
    spectrum(ptbm::ParameterizedTightBindingModel, k::AbstractVector{<:Real})

Evaluate the spectrum, i.e., energies, of the tight-binding model `ptbm` at a single
momentum `k`, across all the bands of `ptbm`.
"""
function spectrum(
    ptbm::ParameterizedTightBindingModel{D},
    k::AbstractVector{<:Real}
) where D
    length(k) == D || error(lazy"dimension mismatch of momentum ($length(k)) & model ($D)")
    H = Hermitian(ptbm(k))
    es = eigvals!(H)
    return es
end