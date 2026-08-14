const HERMITIAN_COLOR = :green
const ANTIHERMITIAN_COLOR = :red

# ---------------------------------------------------------------------------------------- #

function Base.show(io::IO, ::MIME"text/plain", ho::HoppingOrbit)
    # before getting started, determine maximum length of δᵢ entries, to align:
    aligns = map(enumerate(ho.orbit)) do (i, δᵢ)
        s = sprint(print, Crystalline.subscriptify(string(i)), δᵢ)
        textwidth(s)
    end
    max_align = maximum(aligns)
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
        print(io, ": ", " "^(max_align - aligns[i]), "[")
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
        i ≠ length(ho.orbit) && println(io)
    end
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
    length(δs) == 1 && iszero(δs[1]) && return # don't print zero vector (cf. 𝕖(0) = 1)
    if !isnothing(pretext)
        printstyled(io, pretext; stylekws...)
    end
    for (i, δ) in enumerate(δs)
        printstyled(io, "δ", Crystalline.subscriptify(string(i)), "="; stylekws...)
        rev_idx = findfirst(δ′ -> isapprox(-δ, δ′, nothing, false), @view δs[1:i-1])
        if isnothing(rev_idx)
            printstyled(io, replace(string(δ), ", " => ","); stylekws...)
        else
            printstyled(
                io,
                "-δ",
                Crystalline.subscriptify(string(something(rev_idx)));
                stylekws...,
            )
        end
        i == length(δs) || printstyled(io, ", "; stylekws...)
    end
end

# ---------------------------------------------------------------------------------------- #

function _summary_like(io::IO, tbm::TightBindingModel{D, S}, typename::String) where {D, S}
    N = tbm.N
    print(io, length(tbm), "-term ", N, "×", N, " ", typename, "{", D, "}")
    print(io, " (", lowercase(string(S)), ")")
    (N == 0 || length(tbm) == 0) && return
    brs = first(tbm).brs
    print(io, " over ")
    join(io, brs, "⊕")
    return
end
function _summary_like(io::IO, ctbm::CompositeTightBindingModel{D}, typename::String) where D
    tbm_h = ctbm.tbm_h
    tbm_a = ctbm.tbm_a
    N = tbm_h.N
    print(io, "(")
    printstyled(io, length(tbm_h); color=HERMITIAN_COLOR)
    print(io, "+")
    printstyled(io, length(tbm_a); color=ANTIHERMITIAN_COLOR)
    print(io, ")-term ", N, "×", N, " ", typename, "{", D, "}")
    (N == 0 || (length(tbm_h) == 0 && length(tbm_a) == 0)) && return
    brs = first(tbm_h).brs
    print(io, " over ")
    join(io, brs, "⊕")
    return    
end
function Base.summary(io::IO, atbm::AbstractTightBindingModel)
    _summary_like(io, atbm, String(nameof(typeof(atbm))))
end

# This method is factored out of `show(…, ::TightBindingModel)` to allow use also for
# `show(…, ::CompositeTightBindingModel)` without code duplication
function _show_textplain(
    io::IO,
    tbm::TightBindingModel{D};
    add_to_term_counter::Int=0,
    header::Union{Nothing, String}=nothing,
    header_kws=(;),
    trim_line_color=:light_black
) where D
    # NB: Be aware that the first write to `io` in this function is _always_ a newline
    N = tbm.N
    ioc = IOContext(io, :displaysize => displaysize(io) .- (0, 5))
    for (i, tbt) in enumerate(tbm)
        printstyled(io, "\n┌─"; color=trim_line_color)
        if !isnothing(header) && i == 1
            printstyled(io, " ", something(header); header_kws...)
        end
        println(io)
        printstyled(io, i+add_to_term_counter, ". "; bold = true)
        indent = " "^(ndigits(i) + 1)

        s = sprint((io′, x) -> Base.print_array(io′, x), tbt; context = ioc)
        Nˢ = count('\n', s) + 1
        io′ = IOBuffer(s)
        for (i, l) in enumerate(eachline(io′))
            i ≠ 1 && printstyled(io, "│", indent; color=trim_line_color)
            N > 1 && print(io, i == 1 ? '⎡' : i == Nˢ ? '⎣' : '⎢')
            print(io, l)
            N > 1 && print(io, ' ', i == 1 ? '⎤' : i == Nˢ ? '⎦' : '⎥')
            print(io, '\n')
        end
        printstyled(io, "└─ "; color=trim_line_color)
        _print_tightbindingterm_block_summary(io, tbt)
        _print_orbit_elements(io, tbt; color = :light_black, pretext = ":  ")
    end
end
function Base.show(io::IO, ::MIME"text/plain", tbm::TightBindingModel{D}) where {D}
    summary(io, tbm)
    length(tbm) == 0 && return
    print(io, ":")
    return _show_textplain(io, tbm)
end
function Base.show(io::IO, ::MIME"text/plain", ctbm::CompositeTightBindingModel{D}) where {D}
    summary(io, ctbm)
    length(ctbm) == 0 && return
    tbm_h = ctbm.tbm_h
    tbm_a = ctbm.tbm_a
    print(io, ":")
    _show_textplain(
        io, tbm_h;
        header = "Hermitian",
        header_kws = (; color=HERMITIAN_COLOR, italic=true),
        trim_line_color = HERMITIAN_COLOR
    )
    _show_textplain(
        io, tbm_a;
        header = "Anti-Hermitian",
        header_kws = (; color=ANTIHERMITIAN_COLOR, italic=true),
        trim_line_color = ANTIHERMITIAN_COLOR,
        add_to_term_counter = length(tbm_h),
    )
end
function _print_tightbindingterm_block_summary(io::IO, tbt::TightBindingTerm)
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

function Base.summary(io::IO, aptbm::AbstractParameterizedTightBindingModel)
    _summary_like(io, aptbm.tbm, String(nameof(typeof(aptbm))))
end

function Base.show(io::IO, ::MIME"text/plain", aptbm::AbstractParameterizedTightBindingModel)
    summary(io, aptbm)
    length(aptbm.tbm) == 0 && (print(io, " with no amplitudes"); return)
    println(io, " with amplitudes:")
    print(io, " [")
    for (i, c) in enumerate(aptbm.cs)
        if iszero(c)
            printstyled(io, "0"; color = :light_black)
        else
            printstyled(io, round(c, sigdigits=5); color=_coefficient_color(aptbm, i))
        end
        i ≠ length(aptbm.cs) && print(io, ", ")
    end
    print(io, "]")
end
_coefficient_color(::ParameterizedTightBindingModel, _) = :normal
function _coefficient_color(pctbm::ParameterizedCompositeTightBindingModel, i)
    return i ≤ length(pctbm.tbm.tbm_h) ? HERMITIAN_COLOR : ANTIHERMITIAN_COLOR
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
