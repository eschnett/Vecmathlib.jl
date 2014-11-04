module Defs

export @vassert
export vifthen

macro vassert(cond)
end

@inline function vifthen{T}(cond::Bool, x::T, y::T)
    return cond ? x : y
end

end
