function Base.show(io::IO, ::MIME"text/plain", ho::HoppingOrbit)
    # before getting started, determine maximum length of δᵢ entries, and of the associated
    # hopping listings, to align:
    aligns = map(enumerate(ho.orbit)) do (i, δᵢ)
        s = sprint(print, Crystalline.subscriptify(string(i)), δᵢ)
        textwidth(s)
    end
    max_align = maximum(aligns)
    hop_aligns = map(abRs -> textwidth(sprint(_print_hoppings, abRs)), ho.hoppings)
    max_hop_align = maximum(hop_aligns)
    # now print info about each orbit element and its hopping terms
    print(io, typeof(ho), " (")
    printstyled(io, "b"; color = :red)
    print(io, " + ")
    printstyled(io, "R"; color = :blue)
    print(io, " + δ = ")
    printstyled(io, "a"; color = :green)
    println(io, "):")
    for (i, (δᵢ, abRs)) in enumerate(zip(ho.orbit, ho.hoppings))
        print(io, " ")
        printstyled(
            io,
            "δ",
            Crystalline.subscriptify(string(i)),
            " = ",
            δᵢ;
            underline = i == 1,
        )
        print(io)
        print(io, ": ", " "^(max_align - aligns[i]))
        _print_hoppings(io, abRs)
        # if δᵢ is the sign-flipped partner of an earlier orbit element, note the relation:
        # it is what the printed Hamiltonian terms refer to (cf. `_print_orbit_elements`)
        m, negated = canonical_orbit_element(ho.orbit, i)
        if negated
            print(io, " "^(max_hop_align - hop_aligns[i]))
            printstyled(
                io,
                " (= -δ",
                Crystalline.subscriptify(string(m)),
                ")";
                color = :light_black,
            )
        end
        i ≠ length(ho.orbit) && println(io)
    end
end

function _print_hoppings(io::IO, abRs)
    print(io, "[")
    for (j, (a, b, R)) in enumerate(abRs)
        printstyled(io, "("; color = :light_black)
        printstyled(io, b; color = :red)
        print(io, " + ")
        printstyled(io, R; color = :blue)
        print(io, " → ")
        printstyled(io, a; color = :green)
        printstyled(io, ")"; color = :light_black)
        j ≠ length(abRs) && print(io, ", ")
    end
    print(io, "]")
end

# ---------------------------------------------------------------------------------------- #

function Base.summary(io::IO, tbt::TightBindingTerm{D, S}) where {D, S}
    N = last(tbt.axis)
    print(io, N, "×", N, " TightBindingTerm{", D, "}")
    print(io, " (", lowercase(string(S)), ")")
    N == 0 && return
    _print_tightbindingterm_bandreps(io, tbt)
end
function _print_tightbindingterm_bandreps(io::IO, tbt::TightBindingTerm{D}) where {D}
    print(io, " over [")
    for (n, br) in enumerate(tbt.brs)
        printstyled(io, br; color = n ∈ tbt.block_ij ? :blue : :light_black)
        print(io, n == length(tbt.brs) ? "]" : ", ")
    end
end
function Base.show(io::IO, ::MIME"text/plain", tbt::TightBindingTerm{D}) where {D}
    summary(io, tbt)
    println(io, ":")
    ioc = IOContext(io, :displaysize => displaysize(io) .- (0, 3))
    Base.print_array(ioc, tbt)
    _print_orbit_elements(io, tbt; pretext = "\n", color = :light_black, reverse = true)
end

function _print_orbit_elements(
    io::IO,
    tbt::TightBindingTerm;
    pretext = nothing,
    stylekws...,
)
    δs = tbt.block.h_orbit.orbit
    length(δs) == 1 && iszero(δs[1]) && return # don't print zero vector (cf. {0} = 1)
    # only list the canonical representatives of each ±δ pair: their partners are referred
    # to as `-δₘ` in the printed matrix elements, so listing them separately is redundant
    is = filter(i -> !last(canonical_orbit_element(δs, i)), eachindex(δs))
    if !isnothing(pretext)
        printstyled(io, pretext; stylekws...)
    end
    for (n, i) in enumerate(is)
        printstyled(io, "δ", Crystalline.subscriptify(string(i)), "="; stylekws...)
        printstyled(io, replace(string(δs[i]), ", " => ","); stylekws...)
        n == length(is) || printstyled(io, ", "; stylekws...)
    end
end

# ---------------------------------------------------------------------------------------- #

function _summary_like(io::IO, tbm::TightBindingModel{D, S}, spoofname::String) where {D, S}
    N = tbm.N
    print(io, length(tbm), "-term ", N, "×", N, " ", spoofname, "{", D, "}")
    print(io, " (", lowercase(string(S)), ")")
    N == 0 && return
    length(tbm) == 0 && return
    brs = first(tbm).brs
    print(io, " over ")
    join(io, brs, "⊕")
end
Base.summary(io::IO, tbm::TightBindingModel) = _summary_like(io, tbm, "TightBindingModel")

function Base.show(io::IO, ::MIME"text/plain", tbm::TightBindingModel{D}) where {D}
    summary(io, tbm)
    length(tbm) == 0 && return
    print(io, ":")
    N = tbm.N
    ioc = IOContext(io, :displaysize => displaysize(io) .- (0, 5))
    for (i, tbt) in enumerate(tbm)
        printstyled(io, "\n┌─\n"; color = :light_black)
        printstyled(io, i, ". "; bold = true)
        indent = " "^(ndigits(i) + 1)

        s = sprint((io′, x) -> Base.print_array(io′, x), tbt; context = ioc)
        Nˢ = count('\n', s) + 1
        io′ = IOBuffer(s)
        for (i, l) in enumerate(eachline(io′))
            i ≠ 1 && printstyled(io, "│", indent; color = :light_black)
            N > 1 && print(io, i == 1 ? '⎡' : i == Nˢ ? '⎣' : '⎢')
            print(io, l)
            N > 1 && print(io, ' ', i == 1 ? '⎤' : i == Nˢ ? '⎦' : '⎥')
            print(io, '\n')
        end
        printstyled(io, "└─ "; color = :light_black)
        _print_tightbindingterm_block_summary(io, tbt)
        _print_orbit_elements(io, tbt; color = :light_black, pretext = ":  ")
    end
end
function _print_tightbindingterm_block_summary(io::IO, tbt::TightBindingTerm{D}) where {D}
    i, j = tbt.block_ij
    printstyled(io, tbt.brs[i]; color = :blue)
    if i == j
        printstyled(io, " self-term"; color = :light_black)
    else
        print(io, "↔")
        printstyled(io, tbt.brs[j]; color = :blue)
    end
end

# ---------------------------------------------------------------------------------------- #

function Base.summary(io::IO, ptbm::ParameterizedTightBindingModel)
    _summary_like(io, ptbm.tbm, "ParameterizedTightBindingModel")
end

function Base.show(
    io::IO, 
    ::MIME"text/plain", 
    ptbm::ParameterizedTightBindingModel{D}
) where {D}
    summary(io, ptbm)
    length(ptbm.tbm) == 0 && (print(io, " with no amplitudes"); return)
    println(io, " with amplitudes:")
    print(io, " [")
    for (i, c) in enumerate(ptbm.cs)
        if iszero(c)
            printstyled(io, "0"; color = :light_black)
        else
            print(io, round(c, sigdigits=5))
        end
        i ≠ length(ptbm.cs) && print(io, ", ")
    end
    print(io, "]")
end


# ---------------------------------------------------------------------------------------- #

function Base.show(io::IO, cache::TightBindingCache{D, S}) where {D, S}
    N = cache.tbm.N
    Nᵏ = length(cache.ks)
    print(io,
        length(cache.tbm), "-term ", N, "×", N,
        " TightBindingCache{", D, ", …} (", lowercase(string(S)), ") ",
        "over ", Nᵏ, " k-point", Nᵏ == 1 ? "" : "s",
    )
end
