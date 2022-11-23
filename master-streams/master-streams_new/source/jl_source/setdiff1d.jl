using Distributed 
addprocs(8)
using Distributions 
@everywhere using ParallelDataTransfer
@everywhere using SharedArrays
@everywhere using Distributed 
using PythonCall 


# Static number of workers (8 of 'em.)
ws = workers() 

@everywhere function _setdiff1d(
    _set1::SharedVector{Int64}, 
    _set2::SharedVector{Int64}, 
    _set1_indices_to_check::UnitRange{Int64},
    _set1_truefalse::SharedVector{Bool}
    )

    
    @simd for _i in _set1_indices_to_check
        @inbounds if _set1[_i] in _set2 
            @inbounds _set1_truefalse[_i] = true 
        end 
    end   

end 

"""
    _setdiff1d(set1::Vector{Int64}, 
    _set2::Vector{Int64},
    _set1_indices_to_check::Vector{Int64}, 
    _set1_truefalse::SharedVector{Bool})

See above. This is the documentation for the Helper function for setdiff1d, 
the shared array distributed rendition of numpy setdiff1d (but with arguments instead!...)

You should ensure that the cardinality |set1|>=|set2|- I'm not testing otherwise, and it should technically still work, but yeah.
Does not return anything: modifies the shared memory result. 

This particular function makes use of the helper function for distribution.

---
# Arguments 
- `_set1::Vector{Int64}`: the base set which is shared
- `_set2::Vector{Int64}`: the match set which is shared
- `_set1_indices_to_check`: which indices to check for this particular proc (calculated/evaluated)
- `_set1_truefalse::Vector{Bool}`: the result vector
---

"""
function setdiff1d(
    _set1::PyArray{Int64},
    _set2::PyArray{Int64}
)
    _set1 = pyconvert(Array{Int64}, _set1)
    _set2 = pyconvert(Array{Int64}, _set2)

    set1 = SharedVector{Int64}(size(_set1)[1], pids=ws)
    set2 = SharedVector{Int64}(size(_set2)[1], pids=ws)
    set1_truefalse = SharedVector{Bool}(size(set1)[1], pids=ws)

    set1[:] = _set1 
    set2[:] = _set2 
    set1_truefalse .= false 

    # Define the local indices to check
    for i in ws
        local u = @spawnat i localindices(set1)
        u = fetch(u) 
        println("The local indices are... ")
        println(u)
        println("Sending to... ")
        println(i)
        sendto(i, set1_indices_to_check=u)
    end 

    # Distribute and run! 
    futures = Vector{Future}(undef, size(ws)[1])
    for (num,i) in enumerate(ws)
        futures[num] = @spawnat i _setdiff1d(set1,
            set2,
            set1_indices_to_check,
            set1_truefalse
        )
    end 

    # Wait on the futures to all be done...
    for future in futures
        wait(future)
    end 

    # Function is done. Return the result as a regular vector truefalse
    return set1_truefalse[:]

end 
