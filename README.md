Vecmathlib.jl
=============

Vectorizable elemental math functions for Julia

Vecmathlib provides efficient, accurate, tunable, and most importantly
vectorizable math functions such as sqrt, sin, or atan.

# Current State

This repository contains currently only a proof-of-concept
implementation, with code for the `exp` function and a benchmarking
harness. Most ideas are taken from
<https://bitbucket.org/eschnett/vecmathlib>, which is a C++
implementation.

# Benchmark results

**TL;DR:** The exponential function in this library is about twice as
fast as Julia's standard implementation for SIMD-vectorized 64-bit
floating-point operations.

Below are benchmark results from a MacBook Pro 2.7 GHz Intel Core i7
(with AVX instructions). The benchmarks ran on a single core, which
thus likely ran at a higher frequency than 2.7 GHz. The system was
otherwise only lightly used.

Benchmarking paramters were `ni=1000`, i.e. 1000 iterations for an
inner SIMD-parallelized loop, and `nj=1000*1000`, i.e. 1e6 iterations
of this loop. These numbers ensure that all benchmarking data live in
the level 1 data cache. The benchmarking harness performs additional
operations to ensure that these iterations are not optimized away (see
the source code).

All times are in ns (nanoseconds, 1e-9 seconds, smaller is better),
per single amortized function call. That is, this benchmark does not
measure how fast a single call is -- it measures how fast it is to
make many calls in a tight for loop.

Operation | Float32 [ns] | Float64 [ns]
:---------|:-------------|:------------
no-op     | 0.15         | 0.27
add       | 0.17         | 0.31
mul       | 0.17         | 0.31
exp2      | 6.83         | 7.93
vexp2     | 1.08         | 3.39
(yeppp    | ?            | *2.07)

"no-op" performs no operation and measures the overhead of the
benchmarking overhead. "add" and "mul" perform a floating-point
addition and multiplication, respectively. "exp2" is the standard
Julia `exp2` function, "vexp2" is the vectorizable implementation
provided by this library.

The "yeppp" value is an estimate for the performance of the Yeppp
library <http://www.yeppp.info/> according to its documentation, which
lists 5.6 cycles per call for this CPU architecture
<http://www.yeppp.info/home/yeppp-performance-numbers/>. The main
difference in implementation seems to be that Yeppp aggressively
unrolls the SIMD loop, something that Julia/LLVM doesn't do here. (Is
there an `@unroll` macro for Julia?)

As a sanity check, we can compare these numbers to the theoretical
peak performance of this CPU. With AVX instructions, it should execute
4 add and 4 multiply 64-bit operations per cycle, i.e. a single add or
multiply should take 0.09 us plus benchmarking overhead. This is
approximately what we see here. In fact, we measure even a higher
performance, likely because the benchmark payload can be executed in
parallel (superscalar) with the benchmarking harness. That is an
unavoidable measurement error, unless we were to add significantly
more complexity.
