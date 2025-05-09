using GLMakie # TODO: Change to Makie when we convert to package extension

@recipe(HoppingOrbitPlot, h, Rs) do Scene
    Attributes(
        origins = Attributes(
            color = :firebrick2,
            label = "Origins (a)",
        ),
        destinations = Attributes(
            color = :royalblue1,
            label = "Destinations (b+R)",
        ),
        markersize = 0.05,
        bonds = Attributes(
            color = :ivory4,
            linewidth = 2.0,
            label = "Bonds",
        ),
        unitcell = Attributes(
            color = :silver,
            linewidth = 2.0,
            label = "Unit cell",
        ),
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

The `Rs` argument should be a _primitive` basis associated with the lattice underlying `h`.
If omitted, a cubic basis is used.
"""
function Makie.plot!(
    p::HoppingOrbitPlot{<:Tuple{HoppingOrbit{D}, <:DirectBasis{D}}}
) where {D}
    h = p.h[]   # TODO: actually do the Observables updates; just so tedious...
    Rs = p.Rs[]

    P, V = Point{D, Float32}, Vec{D, Float32}
    Rm = stack(Rs)
    origins = mapreduce(vcat, h.hoppings) do hs
        map(hs) do h
            a = Rm*constant(h[1])
            P(a)
        end
    end
    destinations = mapreduce(vcat, h.hoppings) do hs
        map(hs) do h
            b = Rm*constant(h[2])
            R = Rm*constant(h[3])
            P(b + R)
        end
    end

    # plot parallepiped unit cell
    rect = Rect{D, Float32}(P(-0.5), V(1))  # unit cube at origin
    pts = P.(Ref(Rm) .* Makie.GeometryBasics.coordinates(rect))
    if D == 3
        push!(pts, P(NaN))
        unitcell = pts[[1, 3, 4, 2, 1, 5, 6, 2, 9, 6, 8, 7, 5, 9, 8, 4, 9, 7, 3]]
    elseif D == 2
        unitcell = push!(pts, pts[1])
    elseif D == 1
        error("unsupported dimension $D")
    else
        error("unsupported dimension $D")
    end
    lines!(
        p, unitcell; 
        color=p.unitcell.color, linewidth=p.unitcell.linewidth, label=p.unitcell.label
    )

    # plot bonds
    arrows_dir = V.(destinations .- origins)
    arrows!(
        p,
        origins + arrows_dir * 0.125,
        arrows_dir * 0.675;
        color = p.bonds.color, # grayish
        linewidth = D == 3 ? 0.025 : 5.0,
        arrowsize = D == 3 ? V(0.075, 0.075, 0.1) : V(18.75, 25.0),
        label = p.bonds.label,
    )

    # plot atoms
    meshscatter!(
        p,
        unique(destinations);
        markersize = 0.05,
        color = p.destinations.color, label = p.destinations.label
    )
    meshscatter!(
        p,
        unique(origins);
        markersize = 0.05,
        color = p.origins.color, label = p.origins.label
    )

    return p
end

function Makie.convert_arguments(::Type{<:HoppingOrbitPlot}, h::HoppingOrbit{D}) where D
    return (h, _cubic_basis(Val(D)))
end

_cubic_basis(::Val{3}) = crystal(1, 1, 1, π/2, π/2, π/2)
_cubic_basis(::Val{2}) = crystal(1, 1, π/2)
_cubic_basis(::Val{1}) = crystal(1)
_cubic_basis(::Val{D}) where {D} = error("Unsupported dimension: $D")

## --------------------------------------------------------------------------------------- #
# hack overload to set default axis attributes

function Makie.plot(
    h::HoppingOrbit{D},
    Rs::DirectBasis{D} = _cubic_basis(Val(D));
    axis = NamedTuple(),
    figure = NamedTuple(),
    kws...
) where D
    # figure & axis setup
    f = Figure(; figure...)
    ax = if D == 3
        Axis3(f[1, 1]; aspect = :data, viewmode=:fit, axis...)
    else
        Axis(f[1, 1]; aspect = DataAspect(), axis...)
    end
    Makie.hidespines!(ax)
    Makie.hidedecorations!(ax)
    D == 3 && (ax.protrusions[] = 0) # cf. https://github.com/MakieOrg/Makie.jl/issues/2259

    # plot the hopping orbit
    p = hoppingorbitplot!(ax, h, Rs; kws...)

    return Makie.FigureAxisPlot(f, ax, p)
end

Makie.plottype(::HoppingOrbit) = HoppingOrbitPlot
Makie.plottype(::HoppingOrbit{D}, ::DirectBasis{D}) where D = HoppingOrbitPlot
function Makie.args_preferred_axis(
    ::Type{<: HoppingOrbitPlot}, ::HoppingOrbit{D}, ::DirectBasis{D}
) where D
    return D == 3 ? Axis3 : Axis
end
# ---------------------------------------------------------------------------------------- #

# NB: This doesn't give us the nice automatic axis settings that we get from the `plot`
#     overload hack above, but it is the idiomatic way to do it and simpler. If desired, we
#     could overload `Makie.plot` at some point.

Makie.convert_arguments(::Type{<:Plot}, tbb::TightBindingBlock{D}, Rs::DirectBasis{D}) where D = (tbb.h_orbit, Rs)
Makie.convert_arguments(::Type{<:Plot}, tbb::TightBindingBlock{D}) where D = (tbb.h_orbit, _cubic_basis(Val(D)))
Makie.convert_arguments(::Type{<:Plot}, tbt::TightBindingTerm{D}, Rs::DirectBasis{D}) where D = (tbt.block.h_orbit, Rs)
Makie.convert_arguments(::Type{<:Plot}, tbt::TightBindingTerm{D}) where D = (tbt.block.h_orbit, _cubic_basis(Val(D)))