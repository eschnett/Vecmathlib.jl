module FloatProps

export FloatTypes, IntTypes
export floatType, intType
export mkFloat, mkInt, asFloat, asInt
export float_bits
export mantissa_bits, signbit_bits, exponent_bits
export mantissa_mask, exponent_mask, signbit_mask
export exponent_offset, min_exponent, max_exponent



make_mask{T}(nbits::T, nshift::T) =
    convert(T, (convert(T,1) << nbits - convert(T,1)) << nshift)



# Defined in Base:
# eps, nextfloat, precision, prevfloat, realmax, realmin, typemax, typemin
# Inf, Nan

typealias FloatTypes Union(Float16, Float32, Float64)
typealias IntTypes Union(Int16, Int32, Int64)

@inline floatType{T<:FloatTypes}(::Type{T}) = T
@inline intType(::Type{Float16}) = Int16
@inline intType(::Type{Float32}) = Int32
@inline intType(::Type{Float64}) = Int64

@inline mkFloat{T<:FloatTypes}(::Type{T}, x) = convert(T, x)
@inline mkInt{T<:FloatTypes}(::Type{T}, x) = convert(intType(T), x)
@inline asFloat{T<:FloatTypes}(::Type{T}, x) = reinterpret(T, x::intType(T))
@inline asInt{T<:FloatTypes}(::Type{T}, x) = reinterpret(intType(T), x::floatType(T))

@inline float_bits{T<:FloatTypes}(::Type{T}) = mkInt(T, 8*sizeof(T))

@inline mantissa_bits{T<:FloatTypes}(::Type{T}) =
    mkInt(T, precision(convert(T,0)) - 1)
@inline signbit_bits{T<:FloatTypes}(::Type{T}) =
    mkInt(T, 1)
@inline exponent_bits{T<:FloatTypes}(::Type{T}) =
    float_bits(T) - mantissa_bits(T) - signbit_bits(T)

@inline mantissa_mask{T<:FloatTypes}(::Type{T}) =
    make_mask(mantissa_bits(T), mkInt(T,0))
@inline exponent_mask{T<:FloatTypes}(::Type{T}) =
    make_mask(exponent_bits(T), mantissa_bits(T))
@inline signbit_mask{T<:FloatTypes}(::Type{T}) =
    make_mask(signbit_bits(T), exponent_bits(T) + mantissa_bits(T))

@inline exponent_offset{T<:FloatTypes}(::Type{T}) =
    mkInt(T,1) << (exponent_bits(T) - mkInt(T,1)) - mkInt(T,1)
@inline min_exponent{T<:FloatTypes}(::Type{T}) =
    mkInt(T,2) - exponent_offset(T)
@inline max_exponent{T<:FloatTypes}(::Type{T}) =
    mkInt(T,1) << exponent_bits(T) - exponent_offset(T) - mkInt(T,2)

end
