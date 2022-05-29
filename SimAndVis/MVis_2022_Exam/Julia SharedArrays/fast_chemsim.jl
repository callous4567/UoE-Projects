module fast_chemsim 

"""
    oned_threed(__k::Int64, area::Int64, lx::Int64)

Get the (i,j,k) indices corresponding to this contiguous-flattened 1D index, given lx and the corresponding 2D area lx*_ly.
Make sure to specify area/lx as const for performance purposes. Not bounds-safe (will not check for bounds.) 

# Examples
``` 
julia> threed_oned(33, 16, (4,4,4))
(1, 1, 3)
```

"""
@fastmath function oned_threed(__k::Int64, area::Int64, lx::Int64)

    # Get the _k index (zero-based)
    _k = floor(Int64, (__k - 1)/area)

    # Get the _ij index 
    __i = __k - (_k*area)

    # Get the _i index 
    _i = floor(Int64, (__i - 1)/lx)

    # Next the _j index
    _j = floor(Int64, __i - _i*lx)
    
    # Return 
    return _i + 1, _j, _k + 1 

end 

"""
    threed_oned(_i::Int64, _j::Int64, _k::Int64, area::Int64, lx::Int64)

Get the __k contiguous-flattened 1D index corresponding to this (i,j,k), given area and lx. 
This isn't bounds-safe by the way- you might get errors. 

# Examples
```
julia> threed_oned(1,2,3,16,4)
34
```
"""
@fastmath function threed_oned(_i::Int64, _j::Int64, _k::Int64, area::Int64, lx::Int64)
    return (_k-1)*(area) + (_i-1)*lx + _j 
end 


"""
    pbc(_i::Int64, _lx::Int64)

Standard periodic boundary conditions in one dimension.

# Examples
````
julia> pbc(25,20)
5 
````
"""
@fastmath function pbc(_i::Int64, _lx::Int64)

    if _i > (_lx)
        return Int64(1) 
    elseif _i < 1
        return Int64(_lx)
    else 
        return _i 
    end 

end 

"""
    forward_backward_indices(__k::Int64, area::Int64, shape::Tuple{Int64,Int64,Int64})

Return both the forward and backward indices (i+1,j+1,k+1,i-1,j-1,k-1) as tuple for 3D contiguous flattened with PBCs. 

# Examples
```
julia> forward_backward_indices(50, 5, (5,5,10)) 
(30, 46, 75, 45, 49, 25)
```
"""
@fastmath function forward_backward_indices(__k::Int64, area::Int64, shape::Tuple{Int64,Int64,Int64})

    _i,_j,_k = oned_threed(__k, area, shape[1])
    @inbounds _indices = (
        threed_oned(pbc(_i+1,shape[1]), _j, _k, area, shape[1]),
        threed_oned(_i, pbc(_j+1,shape[2]), _k, area, shape[1]),
        threed_oned(_i, _j, pbc(_k+1,shape[3]), area, shape[1]),
        threed_oned(pbc(_i-1,shape[1]), _j, _k, area, shape[1]),
        threed_oned(_i, pbc(_j-1,shape[2]), _k, area, shape[1]),
        threed_oned(_i, _j, pbc(_k-1,shape[3]), area, shape[1])
        )
    return _indices 
end 


"""
    fast_chem(__a::Vector{Float64}, __b::Vector{Float64}, __c::Vector{Float64}, const_1::Float64, const_2::Float64, const_3::Float64, shape::Tuple{Int64,Int64,Int64}, area::Int64, length::Int64)

Iterate forward the three chemical species vectors by one step. Note here that the constants satisfy...

- const_1 = D*dt/dx**2 
- const_2 = q*dt 
- const_3 = p*dt 

The return from this function is a new set of matrices, i.e. 

# Examples 
``` 
julia> fast_chem(a,b,c,c1,c2,c3,shape,length)
__aa,__bb,__cc 
``` 
"""
@fastmath function fast_chem(__a::Vector{Float64}, __b::Vector{Float64}, __c::Vector{Float64}, const_1::Float64, const_2::Float64, const_3::Float64, shape::Tuple{Int64,Int64,Int64}, area::Int64, length::Int64)

    # Pre-allocate 
    __aa, __bb, __cc = Vector{Float64}(undef, length), Vector{Float64}(undef, length), Vector{Float64}(undef, length)

    # Iterate over the array 
    @simd for __i in 1:length

        # Grab the neighbours for the laplacian (indices)
        nbs = forward_backward_indices(__i, area, shape) 

        # Iterate forward the __a vector
        @inbounds __aa[__i] = const_1*(
            __a[nbs[1]] + __a[nbs[2]] + __a[nbs[3]] + __a[nbs[4]] + __a[nbs[5]] + __a[nbs[6]] - 6*__a[__i]
        ) +
        const_2*__a[__i]*(1 - (__a[__i] + __b[__i] + __c[__i])) -
        const_3*__a[__i]*__c[__i] +
        __a[__i] 

        # Iterate forward the __a vector
        @inbounds __bb[__i] = const_1*(
            __b[nbs[1]] + __b[nbs[2]] + __b[nbs[3]] + __b[nbs[4]] + __b[nbs[5]] + __b[nbs[6]] - 6*__b[__i]
        ) +
        const_2*__b[__i]*(1 - (__a[__i] + __b[__i] + __c[__i])) -
        const_3*__a[__i]*__b[__i] +
        __b[__i] 

        # Iterate forward the __a vector
        @inbounds __cc[__i] = const_1*(
            __c[nbs[1]] + __c[nbs[2]] + __c[nbs[3]] + __c[nbs[4]] + __c[nbs[5]] + __c[nbs[6]] - 6*__c[__i]
        ) +
        const_2*__c[__i]*(1 - (__a[__i] + __b[__i] + __c[__i])) -
        const_3*__b[__i]*__c[__i] +
        __c[__i] 
    end 
    
    # Return Vectors
    return __aa, __bb, __cc 

end 

"""
    typefield(a,b,c,indices,length) 

Evaluate the "type" for the chemical vectors a,b,c 
- indices typical (1,2,3,0) const for types (a,b,c,d) dominating with d=1-(a+b+c)

"""
function typefield(__a::Vector{Float64}, __b::Vector{Float64}, __c::Vector{Float64}, indices::Tuple{Int64,Int64,Int64,Int64}, length::Int64)
    __t = Vector{Int64}(undef, length)
    @simd for __i in 1:length 
        @inbounds __t[__i] = indices[argmax((__a[__i], __b[__i], __c[__i], 1-(__a[__i]+__b[__i]+__c[__i])))]
    end 
    return __t 
end 

using IterTools 
"""
    faces(shape::Tuple{Int64,Int64,Int64})

Return the indices for the "visible faces" of a cube when viewed with our usual plotter. 
"""
@fastmath function faces_threed(shape::Tuple{Int64,Int64,Int64})

    # Evaluate the sizes of the tuple vector needed (for partial coords viewing.)
    partial_size = shape[1]*shape[2] + shape[1]*shape[3] + shape[2]*shape[3]
    partial_coords = Vector{Tuple{Int64,Int64,Int64}}(undef, partial_size)

    # Define the indices 
    f1_2 = product(1:shape[1],1:shape[2])
    f3_4 = product(1:shape[1],1:shape[3])
    f5_6 = product(1:shape[2],1:shape[3])
    
    # Iterate over the xy plane 
    for (num,i) in enumerate(product(1:shape[1],1:shape[2]))
        @inbounds partial_coords[num] = (i...,shape[3])
    end 

    # Iterate over the xz plane 
    for (num,i) in enumerate(product(1:shape[1],1:shape[3]))
        @inbounds partial_coords[num + shape[1]*shape[2]] = (i[1],1,i[2])
    end 

    # Iterate over the yz plane 
    for (num,i) in enumerate(product(1:shape[2],1:shape[3]))
        @inbounds partial_coords[num + shape[1]*shape[2] + shape[1]*shape[3]] = (shape[1],i...)
    end 

    """
    f1_2 = bottom_top excess indices (all ij fixed z)  
    f3_4 = sides indices (all iz fixed j) 
    f5_6 = sides indices (all jz fixed i)
    #bot_z = ((d...,1) for d in f1_2)
    #left_x = ((1,d...) for d in f5_6)
    #right_y = ((d[1],shape[2],d[2]) for d in f3_4)
    #fullcoords = (bot_z..., top_z..., left_x..., right_x..., left_y..., right_y...)



    f1_2 = product(1:shape[1],1:shape[2])
    f3_4 = product(1:shape[1],1:shape[3])
    f5_6 = product(1:shape[2],1:shape[3])

    # Grab Indices of all the sides relevant for viewing. 
    top_z = collect(((d...,shape[3]) for d in f1_2))
    right_x = collect(((shape[1],d...) for d in f5_6))
    left_y = collect(((d[1],1,d[2]) for d in f3_4))
    

    # Pre-allocate for the vector of tuples for threed indices 
    partial_fullsize = size(top_z)[1] + size(right_x)[1] + size(left_y)[1]

    # Loop and fill 
    for i in 1:size(top_z)[1]
        partial_coords[i] = top_z[i] 
    end 
    for i in 1:size(right_x)[1]
        partial_coords[i+size(right_x)[1]] = right_x[i] 
    end 
    for i in 1:size(left_y)[1]
        partial_coords[i+size(right_x)[1]+size(left_y)[1]] = left_y[i] 
    end 

    """

    return partial_coords # fullcoords 
end 

"""
    faces_oned(shape::Tuple{Int64,Int64,Int64})

Oned indices in place of the threed indices of faces_threed. Returns both faces_threed AND corresponding oned indices.
"""
@fastmath function faces_oned(shape::Tuple{Int64,Int64,Int64})

    faces_3 = faces_threed(shape) 
    area = shape[1]*shape[2]
    faces_1 = Vector{Int64}(undef, size(faces_3)[1])
    @simd for i in 1:size(faces_3)[1]
        @inbounds faces_1[i] = threed_oned(faces_3[i][1],faces_3[i][2],faces_3[i][3],area,shape[1])
    end 

    return faces_3, faces_1 

end 

using SharedArrays
using Distributed 
@fastmath function Shared_chem(__a::SharedVector{Float64}, __b::SharedVector{Float64}, __c::SharedVector{Float64}, 
    const_1::Float64, const_2::Float64, const_3::Float64, 
    globalshape::Tuple{Int64,Int64,Int64}, globalarea::Int64, locallength::Int64,
    globices::UnitRange{Int64},
    local__aa::Vector{Float64}, local__bb::Vector{Float64}, local__cc::Vector{Float64})

    # Pre-allocate 
    #__aa, __bb, __cc = Vector{Float64}(undef, locallength), Vector{Float64}(undef, locallength), Vector{Float64}(undef, locallength)

    # Iterate over the array 
    @simd for __i in 1:locallength

        # Grab the neighbours for the laplacian (indices)
        nbs = forward_backward_indices(globices[__i], globalarea, globalshape) 

        # Iterate forward the __a vector
        @inbounds local__aa[__i] = const_1*(
            __a[nbs[1]] + __a[nbs[2]] + __a[nbs[3]] + __a[nbs[4]] + __a[nbs[5]] + __a[nbs[6]] - 6*__a[globices[__i]]
        ) +
        const_2*__a[globices[__i]]*(1 - (__a[globices[__i]] + __b[globices[__i]] + __c[globices[__i]])) -
        const_3*__a[globices[__i]]*__c[globices[__i]] +
        __a[globices[__i]] 

        # Iterate forward the __a vector
        @inbounds local__bb[__i] = const_1*(
            __b[nbs[1]] + __b[nbs[2]] + __b[nbs[3]] + __b[nbs[4]] + __b[nbs[5]] + __b[nbs[6]] - 6*__b[globices[__i]]
        ) +
        const_2*__b[globices[__i]]*(1 - (__a[globices[__i]] + __b[globices[__i]] + __c[globices[__i]])) -
        const_3*__a[globices[__i]]*__b[globices[__i]] +
        __b[globices[__i]] 

        # Iterate forward the __a vector
        @inbounds local__cc[__i] = const_1*(
            __c[nbs[1]] + __c[nbs[2]] + __c[nbs[3]] + __c[nbs[4]] + __c[nbs[5]] + __c[nbs[6]] - 6*__c[globices[__i]]
        ) +
        const_2*__c[globices[__i]]*(1 - (__a[globices[__i]] + __b[globices[__i]] + __c[globices[__i]])) -
        const_3*__b[globices[__i]]*__c[globices[__i]] +
        __c[globices[__i]] 
    end 
    
    # Return Vectors
    #return __aa, __bb, __cc 

end 



"""
    shared_typefield(a,b,c,indices,length) 

Evaluate the "type" for the chemical vectors a,b,c 
- indices typical (1,2,3,0) const for types (a,b,c,d) dominating with d=1-(a+b+c)

"""
function shared_typefield(__a::SharedVector{Float64}, __b::SharedVector{Float64}, __c::SharedVector{Float64}, indices::Tuple{Int64,Int64,Int64,Int64}, length::Int64)
    __t = Vector{Int64}(undef, length)
    @simd for __i in 1:length 
        @inbounds __t[__i] = indices[argmax((__a[__i], __b[__i], __c[__i], 1-(__a[__i]+__b[__i]+__c[__i])))]
    end 
    return __t 
end 







end 