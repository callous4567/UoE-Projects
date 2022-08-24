using PythonCall

mas = 1/206264806.2471
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

"""
@fastmath function radecmatch(ra1,dec1,ra2,dec2) # ra1::PyArray{Float64}, dec1::PyArray{Float64}, ra2::PyArray{Float64}, dec2::PyArray{Float64})

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


"""

Get matches between stars using argmin on a distance matrix (basically assuming a 1-1 match exists, guaranteed- close/identical matches 
can cause issues.) 

"""
@fastmath function radecmatch_argmin(ra1::PyArray{Float64}, dec1::PyArray{Float64}, ra2::PyArray{Float64}, dec2::PyArray{Float64})

    ra1 = pyconvert(Array{Float64}, ra1)
    dec1 = pyconvert(Array{Float64}, dec1)
    ra2 = pyconvert(Array{Float64}, ra2)
    dec2 = pyconvert(Array{Float64}, dec2)

    distrix = distance_matrix(ra1, dec1, ra2, dec2)
    matches = ones(Int64, size(ra1)[1])
    for i in 1:size(ra1)[1]
        matches[i] = argmin(distrix[i,:])
    end 

    matches = matches .- 1  
    indices = collect(1:size(matches)[1]) .- 1 
    println("We succeeded in matching! :D Returning the Julia values.")


    return matches, indices 

end 

