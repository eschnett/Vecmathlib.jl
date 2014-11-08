module Log

using Base.Test
using Bitmanip, Defs, FloatProps, Rounding

export vlog, vlog2, vlog10



# Note: The Coefficients were determined by a Mathematica script,
# available at <https://bitbucket.org/eschnett/vecmathlib>. They
# minimize the error in a least-square sense over the interval dx ∈
# [1/sqrt(2); sqrt(2)]. A better alternative would be to minimize the
# maximum of the error.

# One can use a smaller set of coefficients to increase speed at the
# extent of a slightly larger error.

# TODO: Use @evalpoly for this
@inline function vlog2_poly(y2::Float32)
    # Coefficients for Float32,
    # error=7.09807175879142775648452461821e-8
    r = mkFloat(Float32, 0.59723611417135718739797302426)
    r = muladd(r, y2, mkFloat(Float32, 0.961524413175528426101613434))
    r = muladd(r, y2, mkFloat(Float32, 2.88539097665498228703236701))
    return r
end

@inline function vlog2_poly(y2::Float64)
    # Coefficients for Float64,
    # error=1.48294180185938512675770096324e-16
    r = mkFloat(Float64, 0.243683403415639178527756320773)
    r = muladd(r, y2, mkFloat(Float64, 0.26136626803870009948502658))
    r = muladd(r, y2, mkFloat(Float64, 0.320619429891299265439389))
    r = muladd(r, y2, mkFloat(Float64, 0.4121983452028499242926))
    r = muladd(r, y2, mkFloat(Float64, 0.577078017761894161436))
    r = muladd(r, y2, mkFloat(Float64, 0.96179669392233355927))
    r = muladd(r, y2, mkFloat(Float64, 2.8853900817779295236))
    return r
end

const M_SQRT1_2 = 7.071067811865475244008443621048490392848359376884740365883398689953662392310596e-01
const M_SQRT2 = 1.414213562373095048801688724209698078569671875376948073176679737990732478462102e+00
@inline function vlog2{T}(x::T)
    x0 = x
    # Rescale

    e0 = asInt(T,x) & exponent_mask(T) >> mantissa_bits(T)
    e = e0 - exponent_offset(T)
    m = asInt(T,x) & mantissa_mask(T)
    # TODO: handle subnormals

    x = asFloat(T, m | (asInt(T, mkFloat(T,1)) & exponent_mask(T)))
    b = x > M_SQRT2
    x = vifthen(b, mkFloat(T,0.5)*x, x)
    e = vifthen(b, e+mkInt(T,1), e)
    @vassert mkFloat(T,M_SQRT1_2) <= x <= mkFloat(T,M_SQRT2)
    
    y = (x - mkFloat(T,1)) / (x + mkFloat(T,1))
    y2 = y*y
    r = vlog2_poly(y2)
    r *= y
    
    # Undo rescaling
    # TODO: convert to float by
    # - re-interpreting as subnormal, then multiplying with large number
    # - or'ing an exponent, then subtracting the hidden mantissa bit
    r += mkFloat(T, e)
    
    r = vifthen(x0 < mkFloat(T,0) || visnan(x0), mkFloat(T,NaN), r)
    r = vifthen(x0 == mkFloat(T,0), -mkFloat(T,Inf), r)
    r = vifthen(x0 == mkFloat(T,Inf), mkFloat(T,Inf), r)
    
    return r
end

const M_LN2 = 6.931471805599453094172321214581765680755001343602552541206800094933936219696955e-01
@inline function vexp{T}(x::T)
    # TODO: Re-calculate the coefficients for this function to
    # increase accuracy and avoid the multiplication
    return vlog2(x) * mkFloat(T,M_LN2)
end

const M_LN10 = 2.302585092994045684017991454684364207601101488628772976033327900967572609677367e+00
@inline function vlog10{T}(x::T)
    # TODO: Re-calculate the coefficients for this function to
    # increase accuracy and avoid the multiplication
    return vlog2(x) * mkFloat(T,M_LN10)
end



function selftest_type(T::Type)
    zero = convert(T,0)
    min = realmin(T)
    eps = Base.eps(T)
    max = realmax(T)
    vals = T[1.0,
             zero, #TODO nextfloat(zero),
             #TODO prevfloat(min), min, nextfloat(min),
             0.1*eps, eps, 10.0*eps,
             0.1*sqrt(eps), sqrt(eps), 10.0*sqrt(eps),
             1.0-eps, 1.0, 1.0+eps,
             2.0, 10.0, 100.0,
             0.1*sqrt(max), sqrt(max), 10.0*sqrt(max),
             0.1*max, prevfloat(max), max, Inf, NaN]
    for uval in vals
        for val in (uval, -uval)
            local log2_val::T
            try
                log2_val = log2(val)
            catch
                log2_val = convert(T,NaN)
            end
            vlog2_val = vlog2(val)
            @test (isequal(vlog2_val, log2_val) ||
                   isapprox(vlog2_val, log2_val, rtol=eps, atol=eps))
        end
    end
end

function selftest()
    selftest_type(Float32)
    selftest_type(Float64)
end

selftest()

end
