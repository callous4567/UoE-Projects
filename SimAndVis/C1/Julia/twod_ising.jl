# Grab the fast_ising dependency as a module 
cwd = pwd()
include(joinpath(cwd, "fast_ising.jl")) 
using .fast_ising 

# Hard-code parameters for the object (instantiate it)
_lx, _T, _dyn, _binrate, _max_sweeps = Int64(100), Float64(0.1), "g", Int64(250000), Int64(0.7*(10^9))
_temprate = Int64(25*(10^6))
_tempchange = 0.1

# Define _mat (we're going for random instantiation)
using StatsBase 
values = [-1,1]
weights = [1/2,1/2]
_mat = Array{Int64}(undef, (_lx,_lx))
_mat = sample(values, Weights(weights), (_lx,_lx)) 

# Make our object 
_twod_ising = fast_ising.twod_ising(_mat, 
                                    fast_ising.fast_mag(_mat), 
                                    fast_ising.fast_energy(_mat), 
                                    _lx, 
                                    _T, 
                                    _dyn,
                                    _binrate, 
                                    _max_sweeps,
                                    Int64(0),
                                    rand(1:10000000)[1])



# Timing the function (benchmark test)
using BenchmarkTools
function glauberate()
    global _twod_ising 
    fast_ising.fast_glauber(_twod_ising)
end 
@btime glauberate() 



# Animation function 
using Plots 
theme(:dark)
using Printf
anim = Animation(); 
for sweep in collect(0:_twod_ising._max_sweeps)
    global _twod_ising
    if sweep % _twod_ising._binrate == 0
        #save_id = "sweep_" * string(sweep) * ".png"
        temperature = @sprintf("%.2f", _twod_ising._T)
        heatmap(_twod_ising._mat, aspect_ratio=:equal, xlims=[0,_lx], ylims=[0,_lx], dpi=300, legend = :none)
        annotate!((2, 2, text(temperature, :red, :left, :bottom)))
        #savefig(save_id) 
        frame(anim)
    end 
    if sweep % _temprate == 0
        _twod_ising._T += _tempchange 
    end 
    _twod_ising = fast_ising.fast_glauber(_twod_ising)
end 

mp4(anim, "test_anim.mp4", fps=30)  












#using Plots
#scene = heatmap(_mat, aspect_ratio=:equal, xlims=[0,_lx], ylims=[0,_lx], dpi=300, legend = :none) 
#savefig("testplot.png")