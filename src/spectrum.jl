"""
    spectrum(ptbm::ParameterizedTightBindingModel{D, S}, ks; transform = identity)

Evaluate the spectrum, i.e., energies, of the tight-binding model `ptbm` over an iterable
of input momenta `ks`. 

Energies are returned as a matrix, with rows running over momenta and columns over distinct
bands. If `S == HERMITIAN`, the element types (energies) are `Float64`, and otherwise
`ComplexF64`.

## Arguments
- `ptbm`: a [`ParameterizedTightBindingModel`](@ref) to evaluate the spectrum of.
- `ks`: An iterable of momenta. Each such momentum `ks[i]` must evaluate to a real
  `D`-dimensional abstract vector. If `D = 1` (1D models), `ks` can also be any abstract
  vector of real numbers.

## Keyword arguments

- `transform`: a function to apply to the resulting matrix of energies, defaulting to the
  identity function. This can be used to e.g., convert the energies to a different scaling.

## Example

As an example, we evaluating the band structure of graphene.
Below, we first construct and parameterize a tight-binding model for the the (2b|A₁) EBR in
plane group 17, corresponding to the highest-lying orbitals in graphene.
Next, we construct a path along high-symmetry directions of the Brillouin zone using
[Brillouin.jl](https://github.com/thchr/Brillouin.jl), calculate the spectrum across this
path; and finally, plot the band structure using Brillouin and GLMakie (or PlotlyJS):

```julia-repl
julia> using Crystalline, SymmetricTightBinding

julia> brs = calc_bandreps(17, Val(2));

julia> cbr = @composite brs[5]
13-irrep CompositeBandRep{2}:
 (2b|A₁) (2 bands)

julia> ptbm = tb_hamiltonian(cbr, [zeros(Int, dim(cbr))])([0.0, 1.0]);

julia> using Brillouin, GLMakie

julia> kpi = interpolate(irrfbz_path(17, directbasis(17, Val(2))), 100);

julia> plot(kpi, spectrum(ptbm, kpi))
```
"""
function spectrum(
    ptbm::ParameterizedTightBindingModel{D, S},
    ks;
    transform::F = nothing
) where {D, S, F}
    if !(eltype(ks) <: AbstractVector{<:Real})
        error("the elements of `ks` must subtype `AbstractVector{<:Real}`")
    end
    ResultType = S === HERMITIAN ? Float64 : ComplexF64
    Es = Matrix{ResultType}(undef, length(ks), ptbm.tbm.N)
    for (i, k) in enumerate(ks)
        es = spectrum_single_k(ptbm, k; transform)
        @inbounds Es[i, :] .= es
    end
    return Es
end

# if the model is 1D (D==1), we also provide a convenience method that allows passing a
# any abstract vector of real numbers, each an interpreted as independent momentum point
function spectrum(
    ptbm::ParameterizedTightBindingModel{1, S},
    ks::AbstractVector{<:Real},
    transform::F = nothing
) where {S, F}
    ks = [[k] for k in ks]
    return spectrum(ptbm, ks; transform)
end

"""
    spectrum_single_k(ptbm::ParameterizedTightBindingModel, k::AbstractVector{<:Real})

Evaluate the spectrum, i.e., energies, of the tight-binding model `ptbm` at a single
momentum `k`, across all the bands of `ptbm`.
"""
function spectrum_single_k(
    ptbm::ParameterizedTightBindingModel{D, S},
    k::AbstractVector{<:Real};
    transform::F = nothing
) where {D, S, F}
    length(k) == D ||
        error(lazy"dimension mismatch of momentum ($(length(k))) & model ($D)")
    H = ptbm(k)
    es = eigvals!(H)
    return _apply_transform(es, transform)
end

function _apply_transform(Es, transform::F) where F
    if isnothing(transform)
        return Es
    else
        map!(transform, Es, Es)
    end
end
