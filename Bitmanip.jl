module Bitmanip

using FloatProps

export visfinite, visinf, visnan, visnormal, vissubnormal, vsignbit
export vfabs, vflipsign, vcopysign, vneg
export vilogb, vldexp



@inline function visfinite{T}(x::T)
    e = asInt(T,x) & exponent_mask(T)
    return e != exponent_mask(T)
end
@inline function visinf{T}(x::T)
    em = asInt(T,x) & (exponent_mask(T) | mantissa_mask(T))
    return em == exponent_mask(T)
end
@inline function visnan{T}(x::T)
    e = asInt(T,x) & exponent_mask(T)
    m = asInt(T,x) & mantissa_mask(T)
    return e == exponent_mask(T) && m != mkInt(T,0)
end
@inline function visnormal{T}(x::T)
    e = asInt(T,x) & exponent_mask(T)
    return e != exponent_mask(T) && e != mkInt(T,0)
end
@inline function vissubnormal{T}(x::T)
    e = asInt(T,x) & exponent_mask(T)
    m = asInt(T,x) & mantissa_mask(T)
    return e == mkInt(T,0) && m != mkInt(T,0)
end
@inline function vsignbit{T}(x::T)
    s = asInt(T,x) & signbit_mask(T)
    return s != mkInt(T,0)
end



@inline function vfabs{T}(x::T)
    return asFloat(T, asInt(T,x) & ~signbit_mask(T))
end
@inline function vflipsign{T}(x::T, y::T)
    ix = asInt(T,x)
    iysgn = asInt(T,y) & signbit_mask(T)
    return asFloat(T, ix $ iysgn)
end
@inline function vcopysign{T}(x::T, y::T)
    ixmag = asInt(T,x) & ~signbit_mask(T)
    iysgn = asInt(T,y) & signbit_mask(T)
    return asFloat(T, ixmag | iysgn)
end
@inline function vneg{T}(x::T)
    return asFloat(T, asInt(T,x) $ signbit_mask(T))
end



@inline function vilogb{T}(x::T)
    e = asInt(T,x) & exponent_mask(T) >> mantissa_bits(T) - exponent_offset(T)
    return e
end

@inline function vldexp{T}(x::T, i::Integer)
    i::intType(T)
    # Use direct integer manipulation
    # Extract integer as lowest mantissa bits (highest bits still
    # contain offset, exponent, and sign)
    ix = asInt(T,x)
    # Construct scale factor by setting exponent (this shifts out the
    # highest bits)
    scale = asFloat(T, ix << mantissa_bits(T))
    return x * scale
end
@inline function vldexp{T}(x::T, y::T)
    # Use direct integer manipulation
    # Add a large number to shift the integer bits into the rightmost
    # bits of the mantissa. Also already add the exponent offset that
    # we need below.
    offset = mkFloat(T, mkInt(T,1) << mantissa_bits(T) + exponent_offset(T))
    y = y + offset
    # Extract integer as lowest mantissa bits (highest bits still
    # contain offset, exponent, and sign)
    iy = asInt(T,y)
    # Construct scale factor by setting exponent. This shifts out the
    # highest bits, and shifts the lowest mantissa bits into the
    # exponent. We already added the exponent offset above.
    scale = asFloat(T, iy << mantissa_bits(T))
    r = x * scale
    return r
end

end
