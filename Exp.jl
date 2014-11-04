module Exp

using Base.Test
using Bitmanip, Defs, FloatProps, Rounding

export vexp, vexp10, vexp2



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

const M_LOG2E = 1.442695040888963407359924681001892137426645954152985934135449406931109219181187e+00
@inline function vexp{T}(x::T)
    # TODO: Re-calculate the coefficients for this function to
    # increase accuracy and avoid the multiplication
    return vexp2(mkFloat(T,M_LOG2E)*x)
end

const M_LOG10E = 4.342944819032518276511289189166050822943970058036665661144537831658646492088688e-01
@inline function vexp10{T}(x::T)
    # TODO: Re-calculate the coefficients for this function to
    # increase accuracy and avoid the multiplication
    return vexp2(mkFloat(T,M_LOG10E)*x)
end



function selftest_type(T::Type)
    eps = Base.eps(T)
    min = realmin(T)
    max = realmax(T)
    vals = T[0.0, min, nextfloat(min),
             0.1*eps, eps, 10.0*eps,
             0.1*sqrt(eps), sqrt(eps), 10.0*sqrt(eps),
             1.0-eps, 1.0, 1.0+eps,
             2.0, 10.0, 100.0,
             0.1*sqrt(max), sqrt(max), 10.0*sqrt(max),
             0.1*max, prevfloat(max), max, Inf, NaN]
    for uval in vals
        for val in (uval, -uval)
            @test isequal(vexp2(val), exp2(val))
        end
    end
end

function selftest()
    selftest_type(Float32)
    selftest_type(Float64)
end

selftest()

end
