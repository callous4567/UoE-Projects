# Dependencies
using Distributions 

# Grab the fast_xy dependency as a module 
cwd = pwd()
include(joinpath(cwd, "fast_chemsim.jl")) 
using .fast_chemsim 

# Assign constants 
const dx = Float64(1)
const dt = Float64(0.1)
const lx = Int64(512)
const ly = Int64(512)
const lz = Int64(512)
const shape = (lx,ly,lz)
const area = lx*ly 
const length = lx*ly*lz 
const D = Float64(0.5)
const q = Float64(2.5)
const p = Float64(0.5)
current_sweep = Int64(0)
const max_sweeps = Int64(10000)
const binrate = Int64(10)
const c1 = D*dt/(dx^2)
const c2 = q*dt 
const c3 = p*dt 
const indices = (1,2,3,4)

# Generate the initial random matrices for each field (and the typetrix)
a, b, c = rand(Uniform(0, 1/3), shape), rand(Uniform(0, 1/3), shape), rand(Uniform(0, 1/3), shape)
a, b, c = reshape(a, length), reshape(b, length), reshape(c, length)

using BenchmarkTools
@btime fast_chemsim.fast_chem(a,b,c,c1,c2,c3,shape,area,length) 

