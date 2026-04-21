module SymmetricTightBindingMakieExt

## --------------------------------------------------------------------------------------- #

using SymmetricTightBinding
using SymmetricTightBinding: TightBindingTerm, PRUNE_ATOL_DEFAULT
using Crystalline: DirectBasis, crystal, constant, isapproxin
using Makie

## --------------------------------------------------------------------------------------- #

# This method is missing in GeometryBasics (see GeometryBasics.jl PR#277); pirate-patch
# TODO: remove eventually, if above-noted PR is merged and released
function Makie.GeometryBasics.coordinates(rect::Rect{1, T}) where T
    w = rect.widths
    o = rect.origin
    return [Point{1,T}(o[1]), Point{1,T}(o[1]+w[1])]
end

## --------------------------------------------------------------------------------------- #

const MaybeCoefficient = Union{Nothing, <:AbstractVector{<:Number}}
const default_context_attributes = Attributes(;
    include = true, color = :gray75, linecolor = :gray80, linewidth = 1.25, limits = nothing
) # TODO: we define this so we can manually merge: remove this hack once Makie v0.25 is out

@recipe HoppingOrbitPlot (h, Rs, t, offdiag) begin
    origins = Attributes(; color = :firebrick2, label = "Annihilation site (b+R)")
    destinations = Attributes(; color = :royalblue1, label = "Creation sites (a)")
    markersize = 0.05
    bonds = Attributes(; color = :gray37, linewidth = 2.0, label = "Bonds")
    unitcell = Attributes(; color = :gray55, linewidth = 2.0, label = "Unit cell",
                            patchcolor = :gray94)
    context = default_context_attributes
    # TODO: All `Attributes` kwargs are broken (except `context`, which we handle specially)
    #       due to bug in Makie 0.24: currently, caller must give _full_ attributes for any
    #       field that is itself an `Attributes(...)` since no automatic merging of default
    #       attributes with partial caller attributes takes place.
    #       Will apparently be fixed in Makie v0.25+ (cf. ffreyer, Slack).
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
    offdiag = p.offdiag[] # implies that (anti-)hermiticity requires adding reversed hoppings
    
    # TODO: remove manual merging once Makie v0.25 is out (then just `p.context[]` directly)
    context = p.context[]
    if pkgversion(Makie) < v"0.25"
        # Beware that Makie v0.24 merges in the wrong order (i.e., in Makie's `merge(a, b)`,
        # `a` takes precedence over `b`, unlike the usual Julia convention for `merge`)`
        context = merge(context, default_context_attributes)
    end

    P, V = Point{D, Float32}, Vec{D, Float32}
    Rm = stack(Rs)

    # compute creation (`destination`) and annihilation (`origin`) sites for each hopping,
    # accounting for `t` if given
    origins, destinations = if isnothing(t)
        _origins_and_destinations_from_hoppingorbit(h, Rm)
    else
        _origins_and_destinations_from_coefficients(h, t, Rm)
    end

    # For D = 1: lift to 2D for Makie calls; keep D-dimensional originals for context loop
    # Thus, `plot_origins` is either 2D or 3D, while `origins` is D-dimensional; similar for
    # `plot_destinations` and `destinations`.
    plot_origins = D == 1 ? lift_coordinates_to_2D(origins) : origins
    plot_destinations = D == 1 ? lift_coordinates_to_2D(destinations) : destinations

    # plot parallepiped unit cell (with lower left corner at origin)
    rect = Rect{D, Float32}(P(0), V(1)) # unit cube at origin
    pts = P.(Ref(Rm) .* Makie.GeometryBasics.coordinates(rect))
    if D == 3
        push!(pts, P(NaN))
        unitcell = pts[[1, 3, 4, 2, 1, 5, 6, 2, 9, 6, 8, 7, 5, 9, 8, 4, 9, 7, 3]]
    elseif D == 2
        unitcell = pts[[1, 2, 3, 4, 1]]
    elseif D == 1
        # Makie cannot do 1D plotting, so we manually convert the 1D unit cell to a thin
        # 2D unit cell; similarly so, all other coordinates in 1D case are given a y-coord.
        # equal to zero (via `lift_coordinates_to_2D`)
        x_width_1d = if !isnothing(context.limits[])
            (context.limits[] :: Rect{1, Float32}).widths[1]
        else
            max_x = max(maximum(first, origins), maximum(first, destinations), maximum(first, pts))
            min_x = min(minimum(first, origins), minimum(first, destinations), minimum(first, pts))
            max_x - min_x
        end
        y_height = x_width_1d * 0.2f0 # y-height = 20% of x-width
        unitcell = lift_1D_unit_cell_to_2D(pts, y_height) # only vertical lines; NaN-breaks to avoid horizontal lines
        pts = unitcell[[1, 2, 4, 5]] # skip the NaN break so we can draw with `poly!`
    else
        error("unsupported dimension $D")
    end

    if D == 1 || D == 2
        poly!(pts; color=p.unitcell[].patchcolor)
    end
    lines!(
        p,
        unitcell;
        color = p.unitcell[].color,
        linewidth = p.unitcell[].linewidth,
        label = p.unitcell[].label,
    )

    # plot bonds
    if plot_destinations ≠ plot_origins # skip self-energies
        arrows2d_kws = (;
            argmode = :endpoint, # interpret inputs as start/end points, not start/direction
            color = p.bonds[].color, # grayish
            shaftwidth = 3,
            tipwidth = 9,
            tiplength = 10,
            label = p.bonds[].label
        )
        plot_V = D == 1 ? Vec{2, Float32} : V
        arrows_dir = plot_V.(plot_destinations .- plot_origins)
        arrows2d!(
            p,
            plot_origins .+ arrows_dir * .11,
            plot_destinations .- arrows_dir * .09;
            arrows2d_kws...
        )

        # if we're looking at an off-diagonal "block" term, we must also add "reversed"
        # hoppings explicitly, corresponding to the transposed block
        if offdiag
            arrows2d!(
                p,
                plot_destinations .- arrows_dir * 0.11,
                plot_origins .+ arrows_dir .* 0.09;
                arrows2d_kws...
            )
        end
    end

    # plot atoms
    scatter!(
        p,
        unique(plot_destinations);
        marker = :circle,
        markersize = 14,
        color = p.destinations[].color,
        label = p.destinations[].label,
        depth_shift = -0.1,
    )
    scatter!(
        p,
        unique(plot_origins);
        marker = :circle,
        markersize = 14,
        color = p.origins[].color,
        label = p.origins[].label,
        depth_shift = -0.1,
    )

    # set square axis limits, centered around unit cell center
    bbox = if isnothing(context.limits[])
        bbox_coords = Makie.GeometryBasics.coordinates(data_limits(p))
        cntr = sum(Rs) ./ 2
        max_dist = maximum(abs,
            ntuple(d->maximum(v->abs(cntr[d]-getindex(v, d)), bbox_coords), Val(D)))
        pad = maximum(abs, 
            ntuple(d -> splat(-)(extrema(v -> getindex(v, d), filter(!isnan, pts))), Val(D)))
        width = max_dist + pad*0.1
        Rect{D, Float32}(cntr-V(width), V(2*width))
    else
        context.limits[] :: Rect{D, Float32}
    end

    # include symmetry-related sites, that we didn't hop to
    if context.include[]
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
            D ∈ (1, 2) || continue # doesn't look nice for 3D
            _bbox = if D == 1
                # need a 2D bbox to check containment of the 2D-lifted unit cell polygon;
                # full y_height = 2 * half-height, where half-height = abs(unitcell[1][2])
                lift_1D_bbox_to_2D(bbox, 2*abs(unitcell[1][2]))
            else # 2D
                bbox
            end :: Rect{2, Float32}
            _R = D == 2 ? R : Point{2, eltype(R)}(R[1], zero(eltype(R))) # lift to 2D for use in `unitcell′`
            unitcell′ = unitcell .+ Ref(_R)
            any(in(_bbox), unitcell′) || continue # skip if no part is inside bounding box
            lines!(
                p, unitcell′; 
                color=context.linecolor, linewidth=context.linewidth, depth_shift=0
            )
        end
        D == 1 && (rs = lift_coordinates_to_2D(rs))
        scatter!(
            p,
            rs,
            marker = :circle,
            markersize = 11,
            color = context.color,
        )
    end
    if D == 1
        limits!(lift_1D_bbox_to_2D(bbox, 2*abs(unitcell[1][2])))
    else
        limits!(bbox)
    end

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
# 1D to 2D conversion helpers

function lift_1D_unit_cell_to_2D(pts::Vector{Point{1, Float32}}, y_height::Float32)
    x0, x1 = pts[1][1], pts[2][1]
    h = y_height / 2
    # only vertical walls; NaN break avoids drawing horizontal top/bottom lines
    return Point{2, Float32}[
        Point2f(x0, h),  Point2f(x0, -h),
        Point2f(NaN32, NaN32),
        Point2f(x1, -h), Point2f(x1, h),
    ]
end

function lift_1D_bbox_to_2D(bbox::Rect{1, T}, y_height::T) where T
    return Rect{2, T}(bbox.origin[1], -y_height/2, bbox.widths[1], y_height)
end

function lift_coordinates_to_2D(pts::Vector{Point{1, T}}) where T
    return [Point{2, T}(pt[1], zero(T)) for pt in pts]
end

## --------------------------------------------------------------------------------------- #

# Simple case: take no account of a coefficient vector, just plot all hoppings in `h`.
# The hopping `δ` is such that `δ = b+R-a`, so the hopping direction is from `b+R`
# (annihilated) to `a` (created).
function _origins_and_destinations_from_hoppingorbit(
    h::HoppingOrbit{D},
    Rm::AbstractMatrix{<:Real}
) where D
    P = Point{D, Float32}
    destinations = mapreduce(vcat, h.hoppings) do hs
        map(hs) do h
            a = Rm * constant(h[1])
            P(a)
        end
    end
    origins = mapreduce(vcat, h.hoppings) do hs
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
            push!(destinations, a)
            push!(origins, b_plus_R)
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
function _offdiag(tbb::TightBindingBlock{D, S}) where {D, S}
    S === NONHERMITIAN && return false
    return !tbb.diagonal_block
end

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
            plots = [
                S.HoppingOrbitPlot(
                    _orbit(tbt), Rs, _coefficients(tbt), _offdiag(tbt);
                    context = Attributes(; limits = bbox),
                    markersize = 0.1,
                )
            ]
        else
            plots = PlotSpec[]
        end
        ax = if D == 3
            AT(; plots, aspect = :data, viewmode = :fit)
        elseif D == 2 || D == 1
            AT(; plots, aspect = DataAspect())
        else
            error("unsupported dimension D")
        end

        # hiding spines & decorations, in declarative style
        ax.xgridvisible = ax.ygridvisible = false
        ax.xticksvisible = ax.yticksvisible = false
        ax.xticklabelsvisible = ax.yticklabelsvisible = false
        if D == 3
            ax.xspinesvisible = ax.yspinesvisible = ax.zspinesvisible = false
            ax.zgridvisible = ax.zticksvisible = ax.zticklabelsvisible = false
            ax.xlabelvisible = ax.ylabelvisible = ax.zlabelvisible = false
        else # D == 1 or 2
            ax.topspinevisible = ax.bottomspinevisible = false
            ax.leftspinevisible = ax.rightspinevisible = false
            # # disable interactions (don't play nice with x/axislinks)
            # ax.xzoomlock = ax.yzoomlock = ax.xpanlock = ax.ypanlock = true
            # ax.xrectzoom = ax.yrectzoom = false
        end
        ax_bbox = D == 1 ? lift_1D_bbox_to_2D(bbox, bbox.widths[1] * 0.2f0) : bbox
        ax.limits = bbox_to_limits(ax_bbox) # TODO: seems to do nothing?
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
    rect_coords = Makie.GeometryBasics.coordinates(rect)
    sites = P.(Ref(Rm) .* rect_coords)
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
    # TODO: the padding and added width should surely not be the same in all dimensions
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