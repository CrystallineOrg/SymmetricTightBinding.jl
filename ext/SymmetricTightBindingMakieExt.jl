module SymmetricTightBindingMakieExt

## --------------------------------------------------------------------------------------- #

using SymmetricTightBinding
using SymmetricTightBinding: TightBindingTerm, PRUNE_ATOL_DEFAULT
using Crystalline: DirectBasis, crystal, constant, isapproxin
using Makie

## --------------------------------------------------------------------------------------- #

const MaybeCoefficient = Union{Nothing, <:AbstractVector{<:Number}}

@recipe(HoppingOrbitPlot, h, Rs, t, offdiag) do Scene
    Attributes(;
        origins = Attributes(; color = :firebrick2, label = "Origins (a)"),
        destinations = Attributes(; color = :royalblue1, label = "Destinations (b+R)"),
        markersize = 0.05,
        bonds = Attributes(; color = :gray37, linewidth = 2.0, label = "Bonds"),
        unitcell = Attributes(; color = :gray65, linewidth = 2.0, label = "Unit cell",
                                patchcolor = :gray97),
        context = Attributes(; include = true, color = :gray85, 
                               linecolor = :gray90, linewidth= 1.25,
                               limits = nothing),
    )
end

"""
    plot(h::HoppingOrbit{D}, [Rs::DirectBasis{D}])
    plot(tbb::TightBindingBlock{D}, [Rs::DirectBasis{D}])
    plot(tbt::TightBindingTerm{D}, [Rs::DirectBasis{D}])

Visualize the `HoppingOrbit` associated with `h` or a parent structure `tbb` or `tbt` that
embeds a `HoppingOrbit`.
Origin atoms and destination atoms are shown, as well as hoppings between them. The hopping
direction is always from an origin atom to a destination atom, with direction indicated by
an arrow.

The `Rs` argument should be a _primitive_ basis associated with the lattice underlying `h`.
If omitted, a cubic basis is used.
"""
function Makie.plot!(
    p::HoppingOrbitPlot{<:Tuple{
        HoppingOrbit{D},    # h
        DirectBasis{D},     # Rs
        <:MaybeCoefficient, # t
        Bool                # offdiag
    }},
) where {D}
    h = p.h[]   # TODO: actually do the Observables updates; just so tedious...
    Rs = p.Rs[]
    t = p.t[]
    offdiag = p.offdiag[] # implies that (anti-)hermicity requires adding reversed hoppings

    P, V = Point{D, Float32}, Vec{D, Float32}
    Rm = stack(Rs)

    # plot parallepiped unit cell (with lower left corner at origin)
    rect = Rect{D, Float32}(P(0), V(1)) # unit cube at origin
    pts = P.(Ref(Rm) .* Makie.GeometryBasics.coordinates(rect))
    if D == 3
        push!(pts, P(NaN))
        unitcell = pts[[1, 3, 4, 2, 1, 5, 6, 2, 9, 6, 8, 7, 5, 9, 8, 4, 9, 7, 3]]
    elseif D == 2
        unitcell = pts[[1,2,3,4,1]]#push!(pts, pts[1])
    elseif D == 1
        error("unsupported dimension $D")
    else
        error("unsupported dimension $D")
    end
    D == 2 && poly!(pts; color=p.unitcell[].patchcolor)
    lines!(
        p,
        unitcell;
        color = p.unitcell[].color,
        linewidth = p.unitcell[].linewidth,
        label = p.unitcell[].label,
    )

    origins, destinations = if isnothing(t)
        _origins_and_destinations_from_hoppingorbit(h, Rm)
    else
        _origins_and_destinations_from_coefficients(h, t, Rm)
    end


    # plot bonds
    if destinations ≠ origins # skip self-energies
        arrows2d_kws = (;
            argmode = :endpoint, # interpret inputs as start/end points, not start/direction
            color = p.bonds[].color, # grayish
            shaftwidth = 3,
            tipwidth = 9,
            tiplength = 10,
            label = p.bonds[].label
        )
        arrows_dir = V.(destinations .- origins)
        arrows2d!(
            p,
            origins .+ arrows_dir * .11,
            destinations .- arrows_dir * .09;
            arrows2d_kws...
        )
        
        # if we're looking at an off-diagonal "block" term, we must also add "reversed"
        # hoppings explicitly, corresponding to the transposed block
        if offdiag
            arrows2d!(
                p,
                destinations .- arrows_dir * 0.11,
                origins .+ arrows_dir .* 0.09;
                arrows2d_kws...
            )
        end
    end

    # plot atoms
    scatter!(
        p,
        unique(destinations);
        marker = :circle,
        markersize = 14,
        color = p.destinations[].color,
        label = p.destinations[].label,
        depth_shift = -0.1,
    )
    scatter!(
        p,
        unique(origins);
        marker = :circle,
        markersize = 14,
        color = p.origins[].color,
        label = p.origins[].label,
        depth_shift = -0.1,
    )

    # set square axis limits, centered around unit cell center
    bbox = if isnothing(p.context[].limits[])
        bbox_coords = Makie.GeometryBasics.coordinates(data_limits(p))
        cntr = sum(Rs) ./ 2
        max_dist = maximum(abs,
            ntuple(d->maximum(v->abs(cntr[d]-getindex(v, d)), bbox_coords), Val(D)))
        pad = maximum(abs, 
            ntuple(d -> splat(-)(extrema(v -> getindex(v, d), filter(!isnan, pts))), Val(D)))
        width = max_dist + pad*0.1
        Rect{D, Float32}(cntr-V(width), V(2*width))
    else
        p.context[].limits[] :: Rect{D, Float32}
    end

    # include symmetry-related sites, that we didn't hop to
    if p.context[].include[]
        i_lower = Rm \ (bbox.origin)
        i_lower = ntuple(d->floor(Int, i_lower[d])-1, Val(D))
        i_upper = Rm \ (bbox.origin .+ bbox.widths)
        i_upper = ntuple(d->ceil(Int, i_upper[d]), Val(D))
        init_rs = unique!(vcat(origins, destinations))
        rs = Vector{P}()
        for ijk in Iterators.product(ntuple(d->i_lower[d]:i_upper[d], Val(D))...)
            R = Rm * P(ijk)
            for r in init_rs
                r′ = r + R
                r′ in bbox || continue # skip if outside plot limits
                isapproxin(r′, init_rs; atol=1e-5) && continue # skip if previously plotted
                isapproxin(r′, rs; atol=1e-5) && continue      # skip if previously included
                push!(rs, r′)
            end
            # include adjacent unit cell boundaries
            D == 3 && continue # doesn't look nice for 3D
            unitcell′ = unitcell .+ Ref(R)
            any(in(bbox), unitcell′) || continue
            lines!(unitcell′; color=p.context[].linecolor, linewidth=0.5, depth_shift=0)
        end
        scatter!(
            p,
            rs,
            marker = :circle,
            markersize = 11,
            color = p.context[].color,
        )
    end
    limits!(bbox)

    return p
end

function Makie.convert_arguments(::Type{<:HoppingOrbitPlot}, h::HoppingOrbit{D}) where D
    return (h, _cubic_basis(Val(D)), nothing)
end

_cubic_basis(::Val{3}) = crystal(1, 1, 1, π / 2, π / 2, π / 2)
_cubic_basis(::Val{2}) = crystal(1, 1, π / 2)
_cubic_basis(::Val{1}) = crystal(1)
_cubic_basis(::Val{D}) where {D} = error("Unsupported dimension: $D")
## --------------------------------------------------------------------------------------- #

# simple case: take no account of a coefficient vector, just plot all hoppings in `h`
function _origins_and_destinations_from_hoppingorbit(
    h::HoppingOrbit{D},
    Rm::AbstractMatrix{<:Real}
) where D
    P = Point{D, Float32}
    origins = mapreduce(vcat, h.hoppings) do hs
        map(hs) do h
            a = Rm * constant(h[1])
            P(a)
        end
    end
    destinations = mapreduce(vcat, h.hoppings) do hs
        map(hs) do h
            b_plus_R = Rm * (constant(h[2]) + constant(h[3]))
            P(b_plus_R)
        end
    end
    return origins, destinations
end

# complicated case: account for a coefficient vector `t`, which may mean that some hoppings
# in `h` actually do not appear
function _origins_and_destinations_from_coefficients(
    h::HoppingOrbit{D},
    t::AbstractVector{<:Number},
    Rm::AbstractMatrix{<:Real}
) where D
    # reshape `t` to a tensor version `T` with [k, j, i] indices, and `i` denoting `i`th
    # orbit, `j` denoting `j`th hopping in that orbit, and `k` denoting a multi-index into
    # the possible site-irrep-to-site-irrep hoppings associated with `j`
    N = length(t) ÷ 2
    @assert N * 2 == length(t) # storage as [real ..., imag ...]
    I = length(h.hoppings)
    J = length(first(h.hoppings))
    K = N ÷ (I * J)
    @assert I * J * K == N
    Tr = reshape((@view t[1:N]), (K, J, I))
    Ti = reshape((@view t[(N+1):(2N)]), (K, J, I))

    P = Point{D, Float32}
    origins, destinations = Vector{P}(), Vector{P}()
    for (i, hs) in enumerate(h.hoppings)
        # key aim: only add a hopping term if any of its associated coefficients are nonzero
        for (j, h) in enumerate(hs)
            has_hop = (any(x -> abs(x) > PRUNE_ATOL_DEFAULT, @view Tr[:, j, i]) ||
                       any(x -> abs(x) > PRUNE_ATOL_DEFAULT, @view Ti[:, j, i]))
            has_hop || continue
            a = P(Rm * constant(h[1]))
            b_plus_R = P(Rm * (constant(h[2]) + constant(h[3])))
            push!(origins, a)
            push!(destinations, b_plus_R)
        end
    end
    return origins, destinations
end
## --------------------------------------------------------------------------------------- #
# hack overload to set default axis attributes

function Makie.plot(
    h::HoppingOrbit{D},
    Rs::DirectBasis{D} = _cubic_basis(Val(D)),
    t::MaybeCoefficient = nothing,
    offdiag::Bool = false;
    axis = NamedTuple(),
    figure = NamedTuple(),
    kws...,
) where D
    # figure & axis setup
    f = Figure(; figure...)
    ax = if D == 3
        Axis3(f[1, 1]; aspect = :data, viewmode = :fit, axis...)
    else
        Axis(f[1, 1]; aspect = DataAspect(), axis...)
    end
    Makie.hidespines!(ax)
    Makie.hidedecorations!(ax)
    D == 3 && (ax.protrusions[] = 0) # cf. https://github.com/MakieOrg/Makie.jl/issues/2259

    # plot the hopping orbit
    p = hoppingorbitplot!(ax, h, Rs, t, offdiag; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end

Makie.plottype(::HoppingOrbit) = HoppingOrbitPlot
Makie.plottype(::HoppingOrbit{D}, ::DirectBasis{D}) where D = HoppingOrbitPlot
Makie.plottype(::HoppingOrbit{D}, ::DirectBasis{D}, ::MaybeCoefficient) where D = HoppingOrbitPlot
Makie.plottype(::HoppingOrbit{D}, ::DirectBasis{D}, ::MaybeCoefficient, ::Bool) where D = HoppingOrbitPlot
function Makie.args_preferred_axis(
    ::Type{<:HoppingOrbitPlot},
    ::HoppingOrbit{D},
    ::DirectBasis{D},
    ::MaybeCoefficient,
    ::Bool
) where D
    return D == 3 ? Axis3 : Axis
end
# ---------------------------------------------------------------------------------------- #

# NB: This doesn't give us the nice automatic axis settings that we get from the `plot`
#     overload hack above, but it is the idiomatic way to do it and simpler. If desired, we
#     could overload `Makie.plot` at some point.

function Makie.convert_arguments(
    ::Type{<:Plot},
    tbt_or_tbb::Union{TightBindingTerm{D}, TightBindingBlock{D}},
    Rs::DirectBasis{D} = _cubic_basis(Val(D)),
) where D
    (_orbit(tbt_or_tbb), Rs, _coefficients(tbt_or_tbb), _offdiag(tbt_or_tbb))
end

function Makie.plot(
    tbt_or_tbb::Union{TightBindingTerm{D}, TightBindingBlock{D}},
    Rs::DirectBasis{D} = _cubic_basis(Val(D));
    kws...,
) where D
    plot(_orbit(tbt_or_tbb), Rs, _coefficients(tbt_or_tbb), _offdiag(tbt_or_tbb); kws...)
end
_orbit(tbt::TightBindingTerm) = _orbit(tbt.block)
_orbit(tbb::TightBindingBlock) = _orbit(tbb.h_orbit)
_orbit(h::HoppingOrbit) = h
_coefficients(tbt::TightBindingTerm) = _coefficients(tbt.block)
_coefficients(tbb::TightBindingBlock) = tbb.t
_offdiag(tbt::TightBindingTerm) = _offdiag(tbt.block)
_offdiag(tbb::TightBindingBlock) = !tbb.diagonal_block

## --------------------------------------------------------------------------------------- #
# Plotting entire tight-binding models by tiling their hopping terms
import Makie.SpecApi as S

function Makie.convert_arguments(
    ::Type{<:Plot},
    tbm::Union{
            AbstractVector{TightBindingTerm{D}},
            TightBindingModel{D}, 
            ParameterizedTightBindingModel{D}
        },
    Rs::DirectBasis{D} = _cubic_basis(Val(D));
    # kws..., TODO: how to add these with SpecApi?
) where D
    typeof(tbm) === ParameterizedTightBindingModel{D} && (tbm = tbm.tbm)
    Nt = length(tbm)
    n, m = layout_in_grid(Nt)

    bbox = layout_limits(tbm, Rs)
    AT = D == 3 ? S.Axis3 : S.Axis
    axs = D == 3 ? Matrix{typeof(AT())}(undef, n, m) : Matrix{typeof(AT())}(undef, n, m)
    for idx in LinearIndices(axs)
        if idx ≤ length(tbm)
            tbt = tbm[idx]
            plots = [S.HoppingOrbitPlot(_orbit(tbt), Rs, _coefficients(tbt), _offdiag(tbt);
                                        context = (; limits = bbox))]
        else
            plots = PlotSpec[]
        end
        ax = if D == 3
            AT(; plots, aspect = :data, viewmode = :fit)
        else # D == 2
            AT(; plots, aspect = DataAspect())
        end

        # hiding spines & decorations, in declarative style
        ax.xgridvisible = ax.ygridvisible = false
        ax.xticksvisible = ax.yticksvisible = false
        ax.xticklabelsvisible = ax.yticklabelsvisible = false
        if D == 3
            ax.xspinesvisible = ax.yspinesvisible = ax.zspinesvisible = false
            ax.zgridvisible = ax.zticksvisible = ax.zticklabelsvisible = false
            ax.xlabelvisible = ax.ylabelvisible = ax.zlabelvisible = false
        else # D == 2
            ax.topspinevisible = ax.bottomspinevisible = false
            ax.leftspinevisible = ax.rightspinevisible = false
            # # disable interactions (don't play nice with x/axislinks)
            # ax.xzoomlock = ax.yzoomlock = ax.xpanlock = ax.ypanlock = true
            # ax.xrectzoom = ax.yrectzoom = false
        end
        ax.limits = bbox_to_limits(bbox)
        idx ≤ length(tbm) && (ax.title = "Term $idx"; ax.titlefont = :regular; ax.titlegap = 2)
        axs[idx] = ax
    end

    axs = permutedims(axs) # swap from column-major to row-major order for display
    layout = S.GridLayout(axs)
    return layout
end

# find integers `n, m` such that `n * m ≥ p` while `n` and `m` are as close as possible and
# as small as possible, leaving as little "slack" (empty space) as possible
function layout_in_grid(p::Int)
    s = sqrt(p)
    si = round(Int, s)
    si*si == p && return (si, si) # perfect square
    Aₘᵢₙ = typemax(Int)
    Dₘᵢₙ = typemax(Int)
    n′ = m′ = si+1
    for m in 1:ceil(Int, s)
        n = div(p, m, RoundUp)
        A = m * n
        D = abs(m - n)
        if A ≤ Aₘᵢₙ
            if D < Dₘᵢₙ || (D == Dₘᵢₙ && m > m′) # prefer larger m
                n′, m′ = n, m
                Aₘᵢₙ = A
                Dₘᵢₙ = D
            end
        elseif A > Aₘᵢₙ && (D < Dₘᵢₙ && D ≤ div(p, 3, RoundUp))
            n′, m′ = n, m
            Aₘᵢₙ = A
            Dₘᵢₙ = D
        end
    end
    return (n′, m′)
end

function layout_limits(hs::AbstractVector{HoppingOrbit{D}}, Rs::DirectBasis{D}) where D
    P, V = Point{D, Float32}, Vec{D, Float32}
    Rm = stack(Rs)

    # add points from unit cell
    rect = Rect{D, Float32}(P(0), V(1)) # unit cube at origin
    sites = P.(Ref(Rm) .* Makie.GeometryBasics.coordinates(rect))
    filter!(!isnan, sites)
    pts = @view sites[1:end] # for later referencing only the unit cell

    # add points from hoppings
    for h in hs
        for _hs in h.hoppings
            for _h in _hs
                a = P(Rm * constant(_h[1]))                            # origin
                b_plus_R = P(Rm * (constant(_h[2]) + constant(_h[3]))) # destination
                push!(sites, a, b_plus_R)
            end
        end
    end
    unique!(sites)

    # find bounding box
    lower_bound = ntuple(d -> minimum(v -> getindex(v, d), sites), Val(D))
    upper_bound = ntuple(d -> maximum(v -> getindex(v, d), sites), Val(D))
    data_bbox = Rect{D, Float32}(P(lower_bound), V(upper_bound .- lower_bound))
    bbox_coords = Makie.GeometryBasics.coordinates(data_bbox)
    cntr = sum(Rs) ./ 2
    max_dist = maximum(abs,
        ntuple(d->maximum(v->abs(cntr[d]-getindex(v, d)), bbox_coords), Val(D)))
    pad = maximum(abs, 
        ntuple(d -> splat(-)(extrema(v -> getindex(v, d), filter(!isnan, pts))), Val(D)))
    width = max_dist + pad*0.1
    return Rect{D, Float32}(cntr-V(width), V(2*width))
end
function layout_limits(tbm::TightBindingModel{D}, Rs::DirectBasis{D}) where D
    layout_limits([tbt.block.h_orbit for tbt in tbm], Rs)
end
function layout_limits(tbts::AbstractVector{TightBindingTerm{D}}, Rs::DirectBasis{D}) where D
    layout_limits([tbt.block.h_orbit for tbt in tbts], Rs)
end

function bbox_to_limits(bbox::Rect{D}) where D
    lower = bbox.origin
    upper = bbox.origin .+ bbox.widths
    limits = ntuple(Val(2D)) do i # convert to (lower[1], upper[1], lower[2], upper[2], …)
        i′ = (i-1) ÷ 2 + 1
        if mod(i, 2) == 0
            upper[i′]
        else
            lower[i′]
        end
    end
    return limits
end

## --------------------------------------------------------------------------------------- #

end # module SymmetricTightBindingMakieExt