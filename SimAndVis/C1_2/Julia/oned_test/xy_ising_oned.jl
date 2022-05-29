# Grab the fast_xy_oned dependency as a module 
cwd = joinpath(pwd(), "oned_test")
include(joinpath(cwd, "fast_xy_oned.jl")) 
using .fast_xy_oned 
using BenchmarkTools

# Hard-code parameters for the object (instantiate it)
const _lx = Int64(256)
const _T = Float64(0.05)
_binrate = Int64(32768)
_current_bincount = Int64(0)
_pics_per_bin = Int64(90)
_binmultiplier = Int64(2) 
const _max_sweeps = Int64(25*(10^14))

# For temp variation
_temprate = Int64(25*(10^6))
_tempchange = Float64(0.1)


# Define _vec (we're going for random instantiation in 2D to get the angles, then flattening
using Distributions
_vec = rand(Normal(0,1), (_lx, _lx, 2))
for j in 1:_lx 
    for i in 1:_lx 
        _vec[i,j,:] /= sqrt(sum(_vec[i,j,:].^2))
    end
end 
_vec = fast_xy_oned._angle(_vec)
_vec = vcat(_vec...)

"""

# Timing the function (benchmark test)

function glauberate()
    global _vec = fast_xy_oned.fast_onedangleglauber(_vec, _lx, _T)
end 
@btime glauberate() 

"""

# Plotting 
using Printf
using Plots
theme(:dark)
function plot()

    save_id = joinpath(cwd, "BENCHMARK_sweep_" * string(0) * ".png")

    temperature = @sprintf("%.2f", _T)

    heatmap(reshape(_vec, (_lx,_lx)), aspect_ratio=:equal, xlims=[0,_lx], ylims=[0,_lx], dpi=300, clim=(Float64(0), Float64(2*pi)),  c= :thermal)
    
    annotate!((2, 2, text(temperature, :red, :left, :bottom)))
    
    savefig(save_id) 

end 


# Animation function 
for sweep in 0:_max_sweeps

    if sweep % _binrate == 0

        plot()
        
        global _current_bincount += 1 

        if _current_bincount >= _pics_per_bin
            global _binrate *= _binmultiplier
            global _current_bincount = 0 
        end 

    end 

    global _vec = fast_xy_oned.fast_onedangleglauber(_vec, _lx, _T) 

end 
