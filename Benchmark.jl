module Benchmark

using Exp

add1(x) = x + convert(typeof(x), 0.1234)
mul1(x) = x * convert(typeof(x), 0.1234)



function kernel_identity{T}(ni::Int, nj::Int, x::Array{T,1}, y::Array{T,1})
    for j in 1:nj
        @simd for i in 1:ni
            @inbounds y[i] += identity(x[i])
        end
    end
end

function kernel_add{T}(ni::Int, nj::Int, x::Array{T,1}, y::Array{T,1})
    for j in 1:nj
        @simd for i in 1:ni
            @inbounds y[i] += add1(x[i])
        end
    end
end

function kernel_mul{T}(ni::Int, nj::Int, x::Array{T,1}, y::Array{T,1})
    for j in 1:nj
        @simd for i in 1:ni
            @inbounds y[i] += mul1(x[i])
        end
    end
end

function kernel_exp2{T}(ni::Int, nj::Int, x::Array{T,1}, y::Array{T,1})
    for j in 1:nj
        @simd for i in 1:ni
            @inbounds y[i] += exp2(x[i])
        end
    end
end

function kernel_vexp2{T}(ni::Int, nj::Int, x::Array{T,1}, y::Array{T,1})
    for j in 1:nj
        @simd for i in 1:ni
            @inbounds y[i] += vexp2(x[i])
        end
    end
end



function benchmark{T}(func, name, ni::Int, nj::Int, x0::T, dx::T)
    print("Benchmarking $name ($(T)):\n")
    x = Array(T, ni)
    for i in 1:ni
        x[i] = x0 + i*dx
    end
    y = similar(x)
    func(ni, 1, x, y)
    @time func(ni, nj, x, y)
    s = convert(T, 0.0)
    for i in 1:ni
        s += y[i]
    end
    return s
end



function main()
    ni = 1000
    nj = 1000*1000

    x0 = 1.0
    dx = 1.0e-9
    benchmark(kernel_identity, "identity", ni,nj,x0,dx)
    benchmark(kernel_add, "add", ni,nj,x0,dx)
    benchmark(kernel_mul, "mul", ni,nj,x0,dx)
    benchmark(kernel_exp2, "exp2", ni,nj,x0,dx)
    benchmark(kernel_vexp2, "vexp2", ni,nj,x0,dx)

    y0 = 1.0f0
    dy = 1.0f-9
    benchmark(kernel_identity, "identity", ni,nj,y0,dy)
    benchmark(kernel_add, "add", ni,nj,y0,dy)
    benchmark(kernel_mul, "mul", ni,nj,y0,dy)
    benchmark(kernel_exp2, "exp2", ni,nj,y0,dy)
    benchmark(kernel_vexp2, "vexp2", ni,nj,y0,dy)
end

# code_native(kernel_identity, (Int, Int, Array{Float64,1}, Array{Float64,1}))
# code_native(kernel_add, (Int, Int, Array{Float64,1}, Array{Float64,1}))
# code_native(kernel_mul, (Int, Int, Array{Float64,1}, Array{Float64,1}))
# code_native(kernel_exp2, (Int, Int, Array{Float64,1}, Array{Float64,1}))
# code_native(kernel_vexp2, (Int, Int, Array{Float64,1}, Array{Float64,1}))

main()

end
