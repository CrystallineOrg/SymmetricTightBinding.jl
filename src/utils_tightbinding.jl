"""
    split_complex(t::Vector{<:Number}) -> Matrix{Real}

Consider `αt` where `α ∈ ℂ` and `t ∈ ℂⁿ` and build from `t` a matrix representation
`T` that allows access to the real and imaginary parts of the product `αt` without using
complex numbers by splitting α into a real 2-vector of its real and imaginary parts.

In particular, let ``α = αᴿ + iαᴵ`` and ``t = tᴿ + itᴵ`` with `αᴿ, αᴵ ∈ ℝ` and
``tᴿ, tᴵ ∈ ℝⁿ``, then ``αt`` can be rewritten as

```math
αt = (αᴿ + iαᴵ)(tᴿ + itᴵ)
   = (αᴿtᴿ - αᴵtᴵ) + i(αᴿtᴵ + αᴵtᴿ)
   = [tᴿ, tᴵ]ᵀ [αᴿ, αᴵ] + i [tᴵ, tᴿ]ᵀ [αᴿ, αᴵ]
```

Then, defining `T = [tᴿ -tᴵ; tᴵ tᴿ]`, the above product can then be reexpressed as:
``Re(αt) = αᴿtᴿ - αᴵtᴵ =`` `(T * [αᴿ; αᴵ])[1:n]` and ``Im(αt) = αᴿtᴵ + αᴵtᴿ =``
`(T * [αᴿ; αᴵ])[n+1:2n]`.
I.e., the "upper half" of the product `T * [real(α), imag(α)]` is `real(α * t)` and the 
"lower half" is `imag(αt)`.

This functionality is used to avoid complex numbers in amplitude basis coefficients, which
simplifies the application of time-reversal symmetry and hermiticity.

## Examples

```julia
julia> using TETB: split_complex

julia> t = [im,0]
2-element Vector{Complex{Int64}}:
 0 + 1im
 0 + 0im

julia> T = TETB.split_complex(t)
4×2 Matrix{Int64}:
 0  -1
 0   0
 1   0
 0   0

julia> α = 0.5+0.2im; αv = [real(α), imag(α)];

julia> (T * αv)[1:2] == real(α*t) && (T * αv)[3:4] == imag(α*t)
```

```julia
julia> t = [1,im]
2-element Vector{Complex{Int64}}:
 1 + 0im
 0 + 1im

julia> TETB.split_complex(t)
4×2 Matrix{Int64}:
 1   0
 0  -1
 0   1
 1   0
```
"""
function split_complex(t::AbstractVector{<:Number})
    re_t, im_t = reim(t)
    return [re_t -im_t; im_t re_t] # == [real(t) real(im * t); imag(t) imag(im * t)]
end

"""
    inversion(::Val{D}) --> SymOperation{D}

Return the inversion operation in dimension `D`.
"""
inversion(::Val{3}) = S"-x,-y,-z"
inversion(::Val{2}) = S"-x,-y"
inversion(::Val{1}) = S"-x"
inversion(::Val) = error("unsupported dimension")
