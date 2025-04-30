
function Base.summary(io::IO, tbt::TightBindingTerm{D}) where D
    N = last(tbt.axis)
    print(io, N, "×", N, " TightBindingTerm{", D, "}")
    N == 0 && return
    _print_tightbindingterm_bandreps(io, tbt)
end
function _print_tightbindingterm_bandreps(io::IO, tbt::TightBindingTerm{D}) where D
    print(io, " over [")
    for (n, br) in enumerate(tbt.brs)
        printstyled(io, br, color=n ∈ tbt.block_ij ? :blue : :light_black)
        print(io, n == length(tbt.brs) ? "]" : ", ")
    end
end
function Base.show(io::IO, ::MIME"text/plain", tbt::TightBindingTerm{D}) where D
    summary(io, tbt)
    println(io, ":")
    ioc = IOContext(io, :displaysize=>displaysize(io) .- (0,3))
    Base.print_array(ioc, tbt)
    _print_orbit_elements(io, tbt; pretext="\n", color=:light_black, reverse=true)
end

function _print_orbit_elements(io::IO, tbt::TightBindingTerm; pretext=nothing, stylekws...)
    δs = tbt.block.h_orbit.orbit
    length(δs) == 1 && iszero(δs[1]) && return # don't print zero vector (cf. 𝕖(0) = 1)
    if !isnothing(pretext)
        printstyled(io, pretext; stylekws...)
    end
    for (i, δ) in enumerate(δs)
        printstyled(io, "δ", Crystalline.subscriptify(string(i)), "="; stylekws...)
        rev_idx = findfirst(δ′ -> isapprox(-δ, δ′, nothing, false), @view δs[1:i-1])
        if isnothing(rev_idx)
            printstyled(io, replace(string(δ), ", "=>","); stylekws...)
        else
            printstyled(io, "-δ", Crystalline.subscriptify(string(something(rev_idx))); stylekws...)
        end
        i == length(δs) || printstyled(io, ", "; stylekws...)
    end
end

function Base.summary(io::IO, tbm::TightBindingModel{D}) where D
    N = tbm.N
    print(io, length(tbm),"-element ", N, "×", N, " TightBindingModel{", D, "}")
    N == 0 && return
    length(tbm) == 0 && return
    brs = first(tbm).brs
    print(io, " over [")
    for (n, br) in enumerate(brs)
        print(io, br, n == length(brs) ? "]" : ", ")
    end
end

function Base.show(io::IO, ::MIME"text/plain", tbm::TightBindingModel{D}) where D
    summary(io, tbm)
    length(tbm) == 0 && return
    print(io, ":")
    N = tbm.N
    ioc = IOContext(io, :displaysize=>displaysize(io) .- (0,5))
    for (i, tbt) in enumerate(tbm)
        printstyled(io, "\n┌─\n"; color=:light_black)
        printstyled(io, i, ". "; bold=true)
        indent = " "^(ndigits(i) + 1)

        s = sprint((io′, x)->Base.print_array(io′, x), tbt; context=ioc)
        Nˢ = count('\n', s)+1
        io′ = IOBuffer(s)
        for (i, l) in enumerate(eachline(io′))
            i ≠ 1 && printstyled(io, "│", indent, color=:light_black)
            N > 1 && print(io, i==1 ? '⎡' : i==Nˢ ? '⎣' : '⎢')
            print(io, l)
            N > 1 && print(io, ' ', i==1 ? '⎤' : i==Nˢ ? '⎦' : '⎥')
            print(io, '\n')
        end
        printstyled(io, "└─ "; color=:light_black)
        _print_tightbindingterm_block_summary(io, tbt)
        _print_orbit_elements(io, tbt; color=:light_black, pretext=":  ")
    end
end
function _print_tightbindingterm_block_summary(io::IO, tbt::TightBindingTerm{D}) where D
    i, j = tbt.block_ij
    printstyled(io, tbt.brs[i]; color=:blue)
    if i == j
        printstyled(io, " self-term"; color=:light_black)
    else
        print(io, "↔")
        printstyled(io, tbt.brs[j]; color=:blue)
    end
end