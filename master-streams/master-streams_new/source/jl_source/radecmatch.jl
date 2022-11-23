using PythonCall

mas = 1/206264806.2471
"""

Construct a distance matrix in mas such that [i,j] is distance from ra1[i] to ra2[j]. 

"""
@fastmath function distance_matrix(ra1::Array{Float64}, dec1::Array{Float64}, ra2::Array{Float64}, dec2::Array{Float64})

    distrix = Array{Float64}(undef, (size(ra1)[1], size(ra2)[1]))
    @simd for i in 1:size(ra1)[1]
        @simd for j in 1:size(ra2)[1] 
            @inbounds dot_product = cos(dec1[i])*cos(ra1[i])*cos(dec2[j])*cos(ra2[j]) +
                                    cos(dec1[i])*sin(ra1[i])*cos(dec2[j])*sin(ra2[j]) +
                                    sin(dec1[i])*sin(dec2[j])
            if dot_product > 1
                distrix[i,j] = 0 
            elseif dot_product < -1 
                distrix[i,j] = pi 
            else 
                distrix[i,j] = acos(dot_product) 
            end 
        end 
    end 

    return distrix./mas

end 

"""

Construct a vector of length(ra1) with the distances to single float ra2. 
Inverted indices to distance_matrix.

"""
@fastmath function distance_matrix_memlim(ra1::Array{Float64}, dec1::Array{Float64}, ra2::Float64, dec2::Float64)

    distvect = Vector{Float64}(undef, size(ra1)[1])

    @simd for i in 1:size(ra1)[1]

        @inbounds dot_product = cos(dec1[i])*cos(ra1[i])*cos(dec2)*cos(ra2) +
                                cos(dec1[i])*sin(ra1[i])*cos(dec2)*sin(ra2) +
                                sin(dec1[i])*sin(dec2)
        if dot_product > 1
            @inbounds distvect[i] = 0 
        elseif dot_product < -1 
            @inbounds distvect[i] = pi 
        else 
            @inbounds distvect[i] = acos(dot_product) 
        end 
    end 

    return distvect./mas

end 


resolution = 60 # in milliarcseconds 
@fastmath function which_delete(distrix::Array{Float64}, matches::Vector{Tuple{Int64,Int64}})

    to_delete_inverse = ones(Bool, size(matches)[1])
    @simd for i in 1:size(matches)[1]
        @inbounds if distrix[matches[i][1],matches[i][2]] >= resolution 
            to_delete_inverse[i] = false
        end 
    end 

    return to_delete_inverse

end 

using Hungarian 

"""

Get the radec distance matrix for crossmatching using the Munkres algorithm.
    - Best method
    - Minimizes the sum of the excess distances (the cost) as a result of a bad matching
    - Equivalent for argmin in the case of very good matches and no duplicate/neighbours 
Memory expensive (n^2 size arrays.)

"""
@fastmath function radecmatch(ra1::PyArray{Float64}, dec1::PyArray{Float64}, ra2::PyArray{Float64}, dec2::PyArray{Float64})

    ra1 = pyconvert(Array{Float64}, ra1)
    dec1 = pyconvert(Array{Float64}, dec1)
    ra2 = pyconvert(Array{Float64}, ra2)
    dec2 = pyconvert(Array{Float64}, dec2)
    
    distrix = distance_matrix(ra1, dec1, ra2, dec2)
    col_ind, cost = hungarian(distrix)
    matches = Vector{Tuple{Int64,Int64}}(undef, size(ra1)[1])
    distances = ones(Float64, size(matches)[1])
    truefalse = ones(Bool, size(ra1)[1])

    
    @simd for i in 1:size(ra1)[1]
        @inbounds matches[i] = (i, col_ind[i])
        @inbounds matches[i] = matches[i] .- 1 
        @inbounds distances[i] = distrix[i,col_ind[i]]
        @inbounds if distances[i] > resolution
            truefalse[i] = false 
        end 
    end

    println("We succeeded in matching! :D Returning the Julia values.")

    return matches, truefalse

end 

"""

Get matches between stars to within 0.001 mas (assuming that the stars are basically the same ra/dec)
Quick 'n easy but not perfect for dense fields. 

"""
@fastmath function radecmatch_minsep(ra1::PyArray{Float64}, dec1::PyArray{Float64}, ra2::PyArray{Float64}, dec2::PyArray{Float64})

    ra1 = pyconvert(Array{Float64}, ra1)
    dec1 = pyconvert(Array{Float64}, dec1)
    ra2 = pyconvert(Array{Float64}, ra2)
    dec2 = pyconvert(Array{Float64}, dec2)

    matches = ones(Int64, size(ra1)[1])
    found = zeros(Bool, size(ra1)[1])
    @simd for i in 1:size(ra1)[1]
        @simd for j in 1:size(ra2)[1] 
            @inbounds if abs(ra1[i]-ra2[j]) < mas/1000
                @inbounds if abs(dec1[i]-dec2[j]) < mas/1000
                    matches[i] = j 
                    found[i] = true 
                end 
            end 
        end 
    end 

    matches = matches .- 1  

    return matches,found 

end 




@fastmath function radecmatch_argmin(ra1::PyArray{Float64}, dec1::PyArray{Float64}, ra2::PyArray{Float64}, dec2::PyArray{Float64}, argmin_resolution::Int64)

    ra1 = pyconvert(Array{Float64}, ra1)
    dec1 = pyconvert(Array{Float64}, dec1)
    ra2 = pyconvert(Array{Float64}, ra2)
    dec2 = pyconvert(Array{Float64}, dec2)

    distrix = distance_matrix(ra1, dec1, ra2, dec2)
    matches = Vector{Tuple{Int64,Int64}}(undef, size(ra1)[1])
    truefalse = ones(Bool, size(ra1)[1])

    @simd for i in 1:size(ra1)[1]

        # Get argmin (the match)
        col_ind = argmin(distrix[i,:])
        @inbounds matches[i] = (i, col_ind)
        @inbounds matches[i] = matches[i] .- 1 

        # Check if the match exceeds threshold argmin_resolution (in which case reject)
        @inbounds if distrix[i, col_ind] > argmin_resolution
            truefalse[i] = false

        end 
    end 

    println("We succeeded in matching! :D Returning the Julia values.")

    return matches, truefalse 

end 

@fastmath function radecmatch_argmin_memlim(ra1::PyArray{Float64}, dec1::PyArray{Float64}, ra2::PyArray{Float64}, dec2::PyArray{Float64}, argmin_resolution::Int64)

    ra1 = pyconvert(Array{Float64}, ra1)
    dec1 = pyconvert(Array{Float64}, dec1)
    ra2 = pyconvert(Array{Float64}, ra2)
    dec2 = pyconvert(Array{Float64}, dec2)

    matches = Vector{Tuple{Int64,Int64}}(undef, size(ra1)[1])
    distances = zeros(Float64, size(ra1)[1])
    truefalse = ones(Bool, size(ra1)[1])

    @simd for i in 1:size(ra1)[1]

        # Get the distance matrix to this particular ra1 (1-2 inverted here- sorry!)
        distvector = distance_matrix_memlim(ra2, dec2, ra1[i], dec1[i])

        # Get argmin (the match) which is the column index  
        col_ind = argmin(distvector)
        distances[i] = distvector[col_ind]

        # Set them to matches tuple 
        @inbounds matches[i] = (i, col_ind)
        @inbounds matches[i] = matches[i] .- 1 

        # Check if the match exceeds threshold argmin_resolution (in which case reject)
        if distvector[col_ind] > argmin_resolution
            @inbounds truefalse[i] = false

        end 
    end 

    println("We succeeded in matching! :D Returning the Julia values.")

    return matches, distances, truefalse 
end 

