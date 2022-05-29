
"""

# Set dir 
const cwd_local = joinpath(cwd, "test6")
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
using GLMakie, Colors
using GeometryBasics: Rect3D, Point3f0, Vec3f0

# Get indices (threed and oned corresponding) and positions for Makie Meshscatter
const threed_indices, oned_indices = fast_chemsim.faces_oned(shape)
const positions = vec([Point3f0(dd...) for dd in threed_indices])


# Animation function 
for sweep in 0:max_sweeps

    if sweep % binrate == 0

        # Initial typefield for colours 
        typefield = fast_chemsim.typefield(a,b,c,indices,length)
        @inbounds colors = [typefield[d] for d in oned_indices]
        colors[1] = 4 
        
        
        #index = 0 
        #for i in colors
        #    if i == 4
        #        index += 1 
        #    end 
        #end 
        #print(sweep, index, "\n")


        # (Original) by Lazaro Alonso. Figure Production. 
        using GLMakie, Colors, CairoMakie 
        using GeometryBasics: Rect3D, Point3f0, Vec3f0

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

        save(joinpath(cwd_local,"sweep_$current_sweep.png"), fig, px_per_unit=4)
    end 

    if sweep % (50*binrate) == 0 
        GC.gc(true)
    end

    global a,b,c = fast_chemsim.fast_chem(a,b,c,c1,c2,c3,shape,area,length)
    global current_sweep = sweep 
    

end 

"""

"""

fig, ax, _ = volume(x, y, z, vol; colorrange = (1, 4),
    colormap = cs, transparency = true,
    figure = (; resolution = size_in_pix),
    axis = (; type = Axis3, perspectiveness = 0.5, azimuth = 2.19, elevation = 0.57,
        aspect = (1, 1, 1)))
Colorbar(fig[1,2], limits=(1,4), colormap=cs)



# Animation function 
using Plots 
theme(:dark)
using Printf
#anim = Animation(); 
for sweep in 0:max_sweeps

    if sweep % binrate == 0

        save_id = joinpath(cwd_local, "sweep_" * string(sweep) * ".png")

        heatmap(reshape(fast_chemsim.typefield(a,b,c,indices,length), shape)[:,1,:], aspect_ratio=:equal, xlims=[0,lx], ylims=[0,ly], dpi=300, clim=(1,4),  c= :jet)
                
        savefig(save_id) 

        #frame(anim)

    end 


    #if sweep % _temprate == 0
    #    global _T += _tempchange 
    #end  

    global a,b,c = fast_chemsim.fast_chem(a,b,c,c1,c2,c3,shape,area,length)
    global current_sweep = sweep 

end 
"""