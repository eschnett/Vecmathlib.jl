module Benchmark

using Exp

add1(x) = x + convert(typeof(x), 0.1234)
mul1(x) = x * convert(typeof(x), 0.1234)



function kernel_identity{T}(x::Array{T,1}, y::Array{T,1})
    @simd for i in 1:length(y)
        @inbounds y[i] += identity(x[i])
    end
end

function kernel_add{T}(x::Array{T,1}, y::Array{T,1})
    @simd for i in 1:length(y)
        @inbounds y[i] += add1(x[i])
    end
end

function kernel_mul{T}(x::Array{T,1}, y::Array{T,1})
    @simd for i in 1:length(y)
        @inbounds y[i] += mul1(x[i])
    end
end

function kernel_exp2{T}(x::Array{T,1}, y::Array{T,1})
    @simd for i in 1:length(y)
        @inbounds y[i] += exp2(x[i])
    end
end

function kernel_vexp2{T}(x::Array{T,1}, y::Array{T,1})
    @simd for i in 1:length(y)
        @inbounds y[i] += vexp2(x[i])
    end
end



function run_benchmark{T}(func, nj::Int, x::Array{T,1}, y::Array{T,1})
    for j in 1:nj
        func(x, y)
    end
end

function benchmark{T}(func, name, ni::Int, nj::Int, x0::T, dx::T)
    print("Benchmarking $name ($(T)):\n")
    x = Array(T, ni)
    for i in 1:ni
        x[i] = x0 + i*dx
    end
    y = zeros(x)
    run_benchmark(func, 1, x, y)
    @time run_benchmark(func, nj, x, y)
    s = convert(T, 0.0)
    for i in 1:ni
        s += y[i]
    end
    return s
end

function benchmark_type(T::Type)
    ni = 1000
    nj = 1000*1000
    x0 = convert(T, 1.0)
    dx = convert(T, sqrt(eps(T)))
    benchmark(kernel_identity, "identity", ni,nj,x0,dx)
    benchmark(kernel_add, "add", ni,nj,x0,dx)
    benchmark(kernel_mul, "mul", ni,nj,x0,dx)
    benchmark(kernel_exp2, "exp2", ni,nj,x0,dx)
    benchmark(kernel_vexp2, "vexp2", ni,nj,x0,dx)
end


function main()
    benchmark_type(Float32)
    benchmark_type(Float64)
end

# code_native(kernel_identity, (Int, Int, Array{Float64,1}, Array{Float64,1}))
# code_native(kernel_add, (Int, Int, Array{Float64,1}, Array{Float64,1}))
# code_native(kernel_mul, (Int, Int, Array{Float64,1}, Array{Float64,1}))
# code_native(kernel_exp2, (Int, Int, Array{Float64,1}, Array{Float64,1}))
# code_native(kernel_vexp2, (Int, Int, Array{Float64,1}, Array{Float64,1}))

main()

end
