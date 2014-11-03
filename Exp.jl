module Exp

using Base.Test

export vexp, vexp2



macro vassert(cond)
end



# Float32
bits(::Type{Float32}) = int32(32)
floatType(::Type{Float32}) = Float32
intType(::Type{Float32}) = Int32
mkFloat(::Type{Float32}, x) = convert(Float32, x)
mkInt(::Type{Float32}, x) = convert(Int32, x)
asFloat(::Type{Float32}, x) = reinterpret(Float32, x)
asInt(::Type{Float32}, x) = reinterpret(Int32, x)
mantissa_bits(::Type{Float32}) = int32(23)
signbit_bits(::Type{Float32}) = int32(1)
exponent_bits(::Type{Float32}) =
    int32(bits(Float32)) - mantissa_bits(Float32) - signbit_bits(Float32)
mantissa_mask(::Type{Float32}) = int32(1) << mantissa_bits(Float32) - int32(1)
exponent_mask(::Type{Float32}) =
    (int32(1) << exponent_bits(Float32) - int32(1)) << mantissa_bits(Float32)
signbit_mask(::Type{Float32}) =
    (int32(1) << signbit_bits(Float32) - int32(1)) <<
    (mantissa_bits(Float32) + exponent_bits(Float32))
min_exponent(::Type{Float32}) = int32(-125)
max_exponent(::Type{Float32}) =
    int32(1) << exponent_bits(Float32) - exponent_offset(Float32)
exponent_offset(::Type{Float32}) = int32(2) - min_exponent(Float32)

# Float64
bits(::Type{Float64}) = int64(64)
floatType(::Type{Float64}) = Float64
intType(::Type{Float64}) = Int64
mkFloat(::Type{Float64}, x) = convert(Float64, x)
mkInt(::Type{Float64}, x) = convert(Int64, x)
asFloat(::Type{Float64}, x) = reinterpret(Float64, x)
asInt(::Type{Float64}, x) = reinterpret(Int64, x)
mantissa_bits(::Type{Float64}) = int64(52)
signbit_bits(::Type{Float64}) = int64(1)
exponent_bits(::Type{Float64}) =
    int64(bits(Float64)) - mantissa_bits(Float64) - signbit_bits(Float64)
mantissa_mask(::Type{Float64}) = int64(1) << mantissa_bits(Float64) - int64(1)
exponent_mask(::Type{Float64}) =
    (int64(1) << exponent_bits(Float64) - int64(1)) << mantissa_bits(Float64)
signbit_mask(::Type{Float64}) =
    (int64(1) << signbit_bits(Float64) - int64(1)) <<
    (mantissa_bits(Float64) + exponent_bits(Float64))
min_exponent(::Type{Float64}) = int64(-1021)
max_exponent(::Type{Float64}) =
    int64(1) << exponent_bits(Float64) - exponent_offset(Float64)
exponent_offset(::Type{Float64}) = int64(2) - min_exponent(Float64)



# TODO: Use BigFloat to calculate these at run time?
const M_LOG10E = 0.434294481903251827651128918917
const M_LOG2E = 1.44269504088896340735992468100



@inline function vifthen{T}(cond::Bool, x::T, y::T)
    return cond ? x : y
end



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



@inline function vround{T}(x::T)
    # Round by adding, then subtracting again a large number
    offset = mkFloat(T, mkInt(T,1) << mantissa_bits(T))
    tmp = x + offset
    return tmp - offset
end



# Note: The Coefficients were determined by a Mathematica script,
# available at <https://bitbucket.org/eschnett/vecmathlib>. They
# minimize the error in a least-square sense over the interval dx ∈
# [-0.5; +0.5]. A better alternative would be to minimize the maximum
# of the error.

# One can use a smaller set of coefficients to increase speed at the
# extent of a slightly larger error.

# TODO: Use @evalpoly for this
@inline function vexp2_poly(dx::Float32)
    # Coefficients for Float32,
    # error=4.55549108005200277750378992345e-9
    r = mkFloat(Float32, 0.000154653240842602623787395880898)
    r = muladd(r, dx, mkFloat(Float32, 0.00133952915439234389712105060319))
    r = muladd(r, dx, mkFloat(Float32, 0.0096180399118156827664944870552))
    r = muladd(r, dx, mkFloat(Float32, 0.055503406540531310853149866446))
    r = muladd(r, dx, mkFloat(Float32, 0.240226511015459465468737123346))
    r = muladd(r, dx, mkFloat(Float32, 0.69314720007380208630542805293))
    r = muladd(r, dx, mkFloat(Float32, 0.99999999997182023878745628977))
    return r
end

@inline function vexp2_poly(dx::Float64)
    # Coefficients for Float64,
    # error=9.32016781355638010975628074746e-18
    r = mkFloat(Float64, 4.45623165388261696886670014471e-10)
    r = muladd(r, dx, mkFloat(Float64, 7.0733589360775271430968224806e-9))
    r = muladd(r, dx, mkFloat(Float64, 1.01780540270960163558119510246e-7))
    r = muladd(r, dx, mkFloat(Float64, 1.3215437348041505269462510712e-6))
    r = muladd(r, dx, mkFloat(Float64, 0.000015252733849766201174247690629))
    r = muladd(r, dx, mkFloat(Float64, 0.000154035304541242555115696403795))
    r = muladd(r, dx, mkFloat(Float64, 0.00133335581463968601407096905671))
    r = muladd(r, dx, mkFloat(Float64, 0.0096181291075949686712855561931))
    r = muladd(r, dx, mkFloat(Float64, 0.055504108664821672870565883052))
    r = muladd(r, dx, mkFloat(Float64, 0.240226506959101382690753994082))
    r = muladd(r, dx, mkFloat(Float64, 0.69314718055994530864272481773))
    r = muladd(r, dx, mkFloat(Float64, 0.9999999999999999978508676375))
    return r
end

@inline function vexp2{T}(x::T)
    # Here it does not matter how ties are broken
    xi = vround(x)
    dx = x - xi
    @vassert mkFloat(T,-0.5) <= dx <= mkFloat(T,0.5)
    r = vexp2_poly(dx)
    # r = vldexp(r, mkInt(T, xi))
    r = vldexp(r, xi)
    r = vifthen(x < mkFloat(T, min_exponent(T)), mkFloat(T,0.0), r)
    r = vifthen(x > mkFloat(T, max_exponent(T)), mkFloat(T,Inf), r)
    return r
end

@inline function vexp{T}(x::T)
    # TODO: Re-calculate the coefficients for this function to
    # increase accuracy and avoid the multiplication
    return vexp2(mkFloat(T,M_LOG2E)*x)
end

@inline function vexp10{T}(x::T)
    # TODO: Re-calculate the coefficients for this function to
    # increase accuracy and avoid the multiplication
    return vexp2(mkFloat(T,M_LOG10E)*x)
end



function selftests()
    eps32 = eps(Float32)
    min32 = realmin(Float32)
    max32 = realmax(Float32)
    vals32 = Float32[0.0, min32, nextfloat(min32),
                     0.1*eps32, eps32, 10.0*eps32,
                     0.1*sqrt(eps32), sqrt(eps32), 10.0*sqrt(eps32),
                     1.0-eps32, 1.0, 1.0+eps32,
                     2.0, 10.0, 100.0,
                     0.1*sqrt(max32), sqrt(max32), 10.0*sqrt(max32),
                     0.1*max32, prevfloat(max32), max32, Inf, NaN]
    for uval in vals32
        for val in (uval, -uval)
            @test isequal(vexp2(val), exp2(val))
        end
    end
    
    eps64 = eps(Float64)
    min64 = realmin(Float64)
    max64 = realmax(Float64)
    vals64 = Float64[0.0, min64, nextfloat(min64),
                     0.1*eps64, eps64, 10.0*eps64,
                     0.1*sqrt(eps64), sqrt(eps64), 10.0*sqrt(eps64),
                     1.0-eps64, 1.0, 1.0+eps64,
                     2.0, 10.0, 100.0,
                     0.1*sqrt(max64), sqrt(max64), 10.0*sqrt(max64),
                     0.1*max64, prevfloat(max64), max64, Inf, NaN]
    for uval in vals64
        for val in (uval, -uval)
            @test isequal(vexp2(val), exp2(val))
        end
    end
end

selftests()

end
