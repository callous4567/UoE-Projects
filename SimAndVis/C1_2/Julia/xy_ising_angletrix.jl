# Grab the fast_xy dependency as a module 
cwd = pwd()
include(joinpath(cwd, "fast_xy.jl")) 
using .fast_xy 

# Hard-code parameters for the object (instantiate it)
const _lx = Int64(256)
const _T = Float64(0.05)
_binrate = Int64(1)
_current_bincount = Int64(0)
_pics_per_bin = Int64(90)
_binmultiplier = Int64(2) 
const _max_sweeps = Int64(25*(10^14))

# For temp variation
_temprate = Int64(25*(10^6))
_tempchange = Float64(0.1)


# Define _mat (we're going for random instantiation in 2D)
using Distributions
_mat = rand(Normal(0,1), (_lx, _lx, 2))
for j in 1:_lx 
    for i in 1:_lx 
        _mat[i,j,:] /= sqrt(sum(_mat[i,j,:].^2))
    end
end 
_mat = fast_xy._angle(_mat)

# Convert _mat into a vector of vectors
_vec = Vector{Vector{Float64}}(undef, _lx)
for i in 1:_lx
    _vec[i] = _mat[i,:]
end 




# Timing the function (benchmark test)

using BenchmarkTools
function glauberate()
    global _mat = fast_xy.fast_angleglauber(_mat, _lx, _T)
end 
@btime glauberate() 
"""


# Animation function 
using Plots 
theme(:dark)
using Printf
#anim = Animation(); 
for sweep in 0:_max_sweeps

    if sweep % _binrate == 0

        save_id = "sweep_" * string(sweep) * ".png"

        temperature = @sprintf("%.2f", _T)

        heatmap(_mat, aspect_ratio=:equal, xlims=[0,_lx], ylims=[0,_lx], dpi=300, clim=(Float64(0), Float64(2*pi)),  c= :thermal)
        
        annotate!((2, 2, text(temperature, :red, :left, :bottom)))
        
        savefig(save_id) 
        
        global _current_bincount += 1 

        if _current_bincount >= _pics_per_bin
            global _binrate *= _binmultiplier
            global _current_bincount = 0 
        end 

        #frame(anim)

    end 


    #if sweep % _temprate == 0
    #    global _T += _tempchange 
    #end  


    global _mat = fast_xy.fast_angleglauber(_mat, _lx, _T) 

end 

#mp4(anim, "test_anim.mp4", fps=30)  

"""