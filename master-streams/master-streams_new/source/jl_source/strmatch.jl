using PythonCall

"""

Find the indices in list2 corresponding to list1, assuming they exist. Return is sparse.

    # This is for PyArrays only. 
    # We assume the indices *DO EXIST* 

"""
@fastmath function strmatch(list1, list2)

    list1 = pyconvert(Array{String}, list1)
    list2 = pyconvert(Array{String}, list2)
    matches = ones(Int64, size(list1)[1])
    @simd for j in 1:size(list2)[1]
        @simd for i in 1:size(list1)[1] 
            @inbounds if list2[j] == list1[i] 
                matches[i] = j 
            end 
        end 
    end 

    matches = matches .- 1 

    return matches 

end 


"""

Find the indices in list2 corresponding to list1, allowing for non-existence 

"""
@fastmath function strmatch_general(list1, list2)

    list1 = pyconvert(Array{String}, list1)
    list2 = pyconvert(Array{String}, list2)
    matches = Array{Int64}(undef, size(list1)[1])
    found = zeros(Bool, size(list1)[1])
    @simd for j in 1:size(list2)[1]
        @simd for i in 1:size(list1)[1] 
            @inbounds if list2[j] == list1[i] 
                matches[i] = j 
                found[i] = true 
            end 
        end 
    end 

    matches = matches .- 1 

    return matches, found

end 