using Plots

# Define the x-range 
x_range = Array{Float64}(undef, (1,2))
x_range[1],x_range[2] = Float64(0), Float64(8)

# Define the y-range 
y_range = Array{Float64}(undef, (1,2))
y_range[1], y_range[2] = Float64(-4), Float64(2)

# Define the z-range 
z_range = Array{Float64}(undef, (1,2))
z_range[1], z_range[2] = Float64(-4), Float64(2)

# Define the t-range 
num_steps = Int64(1000)
t_range = collect(LinRange(Float64(0), num_steps, num_steps))

# Generate linspaces and collect the arrays (need to collect after instantiating object, or else you just get a generator)
xx = collect(LinRange(x_range[1], x_range[2], Int64(10000)))

# Constants for functions 
k1, w1, d1 = Float64(1), Float64(0.1), Float64(pi/2)
k2, w2, d2 = Float64(1), Float64(0.2), Float64(0) 

# Create a zeros array for if needed
_zeros = fill(Float64(0), size(xx))

# Evaluate y(x,t)
function y_1(_xx::Vector{Float64},t::Float64)

    _yy = sin.(k1*_xx .- w1*t .+ d1)
    return _yy 

end 

# Evaluate z(x,t)
function z_1(_xx::Vector{Float64},t::Float64)
    _zz = sin.(k2*_xx .- w2*t .+ d2)
    return _zz

end 

# Animation function 
using Plots 
anim = Animation(); 
for t in t_range
    _size = size(xx)
    y_columns = Array{Float64}(undef, (_size[1], 3))
    y_columns[:,1], y_columns[:,2], y_columns[:,3] = y_1(xx, t), _zeros, y_1(xx, t)
    z_columns = Array{Float64}(undef, (_size[1], 3))
    z_columns[:,1], z_columns[:,2], z_columns[:,3] = _zeros, z_1(xx, t), z_1(xx, t)
    plot(xx, y_columns, z_columns, legend=false, dpi=300)
    frame(anim)
end 

mp4(anim, "test_anim.mp4", fps=15)
