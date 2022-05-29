# Grab the fast_xy dependency as a module 
cwd = pwd()
include(joinpath(cwd, "fast_xy.jl")) 
using .fast_xy 

# Hard-code parameters for the object (instantiate it)
_lx, _T, _dyn, _binrate, _max_sweeps = Int64(256), Float64(0.05), "g", Int64(2500000), Int64(25*(10^14))
_identifier = rand(1:10000000)[1]
_current_sweep = Int64(0)

# For temp variation
_temprate = Int64(25*(10^6))
_tempchange = 0.1

# Define _mat (we're going for random instantiation in 2D)
using StatsBase 
using Distributions
_mat = rand(Normal(0,1), (_lx, _lx, 2))
for j in 1:_lx 
    for i in 1:_lx 
        _mat[i,j,:] /= sqrt(sum(_mat[i,j,:].^2))
    end
end 
_M, _E = fast_xy.fast_mag(_mat), fast_xy.fast_energy(_mat)



# Timing the function (benchmark test)
using BenchmarkTools
function glauberate()
    fast_xy.fast_glauber(_mat, _lx, _M, _E, _T)
end 
@btime glauberate() 



# Animation function 
using Plots 
theme(:dark)
using Printf
#anim = Animation(); 
for sweep in 0:_max_sweeps

    if sweep % _binrate == 0
        save_id = "sweep_" * string(sweep) * ".png"
        temperature = @sprintf("%.2f", _T)
        angle_mat = fast_xy._angle(_mat)
        heatmap(angle_mat, aspect_ratio=:equal, xlims=[0,_lx], ylims=[0,_lx], dpi=300, c= :thermal, clim=(0, 2*pi))
        annotate!((2, 2, text(temperature, :red, :left, :bottom)))
        savefig(save_id) 
        #frame(anim)

    end 

    """
    if sweep % _temprate == 0
        global _T += _tempchange 
    end  

    """

    global _mat, _M, _E = fast_xy.fast_glauber(_mat, _lx, _M, _E, _T) 

end 

#mp4(anim, "test_anim.mp4", fps=30)  
