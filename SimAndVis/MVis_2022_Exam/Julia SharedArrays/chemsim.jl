using Distributed 
addprocs(10)
using Distributions 
@everywhere using ParallelDataTransfer
@everywhere using SharedArrays
@everywhere using Distributed 

# Static number of workers (arranged for half/half plot/simulate)
ws = workers() 
ws_plot = ws[1:8]
ws_sim = ws[9:10] 

# Grab the fast_xy dependency as a module 
@everywhere cwd = pwd()
@everywhere include(joinpath(cwd, "fast_chemsim.jl")) 
@everywhere using .fast_chemsim 
@everywhere include(joinpath(cwd, "chemplot.jl"))
@everywhere using .chemplot 

# Assign constants 
@everywhere const dx = Float64(1)
@everywhere const dt = Float64(0.1)
@everywhere const lx = Int64(256)
@everywhere const ly = Int64(256)
@everywhere const lz = Int64(256)
@everywhere const shape = (lx,ly,lz)
@everywhere const area = lx*ly 
@everywhere const length = lx*ly*lz 
@everywhere const D = Float64(1)
@everywhere const q = Float64(0.5)
@everywhere const p = Float64(4)
current_sweep = Int64(0)
@everywhere const max_sweeps = Int64(10000)
@everywhere const binrate = Int64(10)
@everywhere const c1 = D*dt/(dx^2)
@everywhere const c2 = q*dt 
@everywhere const c3 = p*dt 
@everywhere const indices = (1,2,3,4)


# Generate the initial random matrices for each field 
a, b, c = rand(Uniform(0, 1/3), shape), rand(Uniform(0, 1/3), shape), rand(Uniform(0, 1/3), shape)
a, b, c = reshape(a, length), reshape(b, length), reshape(c, length)
aa, bb, cc = SharedArray{Float64}(size(a), pids=ws_sim),
    SharedArray{Float64}(size(a), pids=ws_sim),
    SharedArray{Float64}(size(a), pids=ws_sim)
aa[:], bb[:], cc[:] = a,b,c 

# Make sure to pre-define the ranges for simulation indices 
for (num,i) in enumerate(ws_sim) 
    local u = @spawnat i localindices(aa)
    u = fetch(u) 
    sendto(i, local_globices=u)
end 

# Define the local sets in all 
for i in ws_sim 

    @defineat i begin 

        locallength = size(local_globices)[1]
        local__aa, local__bb, local__cc = Vector{Float64}(undef, locallength), Vector{Float64}(undef, locallength), Vector{Float64}(undef, locallength)

    end 

end 

# Function to change the locals 
@everywhere function update_locals(aa,bb,cc)
    fast_chemsim.Shared_chem(aa, bb, cc, 
        c1, c2, c3, 
        shape, area, 
        size(local_globices)[1], local_globices,
        local__aa, local__bb, local__cc)
    # global local__aa, local__bb, local__cc = 
end 

# Function to set the shared array 
@everywhere function set_shared(aa,bb,cc)
    aa[local_globices], 
    bb[local_globices], 
    cc[local_globices] = local__aa, local__bb, local__cc
end 

function iter_step()

    # Store the Futures to wait on
    futures = Vector{Future}(undef, size(ws_sim)[1])

    # Run the shared_chem steps 
    for (num,i) in enumerate(ws_sim) 
        futures[num] = @spawnat i update_locals(aa,bb,cc)
    end 

    # Wait on all of them to finish updating localsy 
    for future in futures 
        wait(future)
    end 

    # All prepared. Set all locals 
    for (num,i) in enumerate(ws_sim) 
        futures[num] = @spawnat i set_shared(aa,bb,cc)
    end 

    # Wait on all of them to finish updating localsy 
    for future in futures 
        wait(future)
    end 
end 

using BenchmarkTools
@btime iter_step() 

# Define a placeholder to keep track of all our spawns 
spawns = Vector{Future}(undef,0)
last_process = 0

# Animation function 
for sweep in 0:max_sweeps

    if sweep % binrate == 0

        # Specify Typefield 
        typefield = fast_chemsim.shared_typefield(aa,bb,cc,indices,length)

        # Push typefield on to new process (cycling through processes in ws_sim)
        if last_process == 8 
            global last_process = 0
        end 
        pushed = @spawnat ws_plot[last_process+1] chemplot.plot(typefield, sweep)
        global last_process += 1 
        push!(spawns, pushed)

        # If the length of spawns exceeds 16, try fetching the first half and clearing them
        if size(spawns)[1] > 30 
            
            for i in 1:15
                fetch(spawns[i])
            end 

            println("Fetched all spawns.")

            global spawns = spawns[9:size(spawns)[1]]

        end 


    end 

    iter_step()
    global current_sweep = sweep 

    if sweep % (10*binrate) == 0 
        @everywhere GC.gc(true)
    end 
    

end 