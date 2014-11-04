module Rounding

using FloatProps

export vround



@inline function vround{T}(x::T)
    # Round by adding, then subtracting again a large number
    offset = mkFloat(T, mkInt(T,1) << mantissa_bits(T))
    tmp = x + offset
    return tmp - offset
end

end
