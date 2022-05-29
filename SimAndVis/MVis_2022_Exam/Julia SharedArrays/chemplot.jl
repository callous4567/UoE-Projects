module chemplot 

# Grab the fast_chemsim dependency as a module 
cwd = pwd()
include(joinpath(cwd, "fast_chemsim.jl")) 
using .fast_chemsim 

# Assign constants 
const dx = Float64(1)
const dt = Float64(0.1)
const lx = Int64(256)
const ly = Int64(256)
const lz = Int64(256)
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

# Set dir 
const cwd_local = joinpath(cwd, "test_multicore_5")
try 
    mkdir(cwd_local) 
catch 
    print("Directory exists \n") 
end 

# Make a colourmap 
using Colors, ColorSchemes
c11, c22, c33, c44 = colorant"red", colorant"green", colorant"blue", colorant"gray"
const cs = ColorScheme([c11,c22,c33,c44])


# Plotting 
using GLMakie, Colors, CairoMakie
using GeometryBasics: Rect3D, Point3f0, Vec3f0

# Get indices (threed and oned corresponding) and positions for Makie Meshscatter
const threed_indices, oned_indices = fast_chemsim.faces_oned(shape)
const positions = vec([Point3f0(dd...) for dd in threed_indices])


function plot(typefield::Vector{Int64}, sweep::Int64)

    @inbounds colors = [typefield[d] for d in oned_indices]
    colors[1] = 4 
    colors[2] = 3 
    colors[3] = 2
    colors[4] = 1 

    # (Original) by Lazaro Alonso. Figure Production. 

    fig, ax, obj = meshscatter(positions; marker = Rect3D(Vec3f0(-8), Vec3f0(8)),
        transparency = false,
        shading = false,
        color = colors,
        figure = (;),
        colormap=cs,
        clim=(0.9,4.1),
        axis = (; type = Axis3, perspectiveness = 0.5, azimuth = -1*pi/4, elevation = pi/4,
                aspect = (1,1,1))
    )
    Colorbar(fig[1,2], limits=(0.9,4.1), colormap=cs)
    GLMakie.resize!

    save(joinpath(cwd_local,"sweep_$sweep.png"), fig, px_per_unit=4)

    if sweep % (50*binrate) == 0 
        GC.gc(true)
    end 

end 

end 