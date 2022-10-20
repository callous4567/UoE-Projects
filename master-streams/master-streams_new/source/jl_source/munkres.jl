using Hungarian 

function solve_hungarian(cost_matrix::Array{Int64,2})

    assignment, cost = hungarian(cost_matrix)

    return assignment 

end 

using FreqTables
using PythonCall 

"""

    compclust_multilabel(clustt1::Vector{Int64}, clustt2::Vector{Int64}, minimum_excess_index::Int64)

Remap clust2 to clust1 with unmapped clusters in clust2 being relabelled sequentially starting from minimum_excess_index. 

# Examples 
``` 
julia> a = [1, 1, 2, 2, 2, 2]
julia> b = [2, 2, 3, 3, 6, 2]
julia> compclust_multilabel(a,b,100)
[1, 1, 2, 2, 
```

"""
@fastmath function compclust_multilabel(clustt1::PyArray{Int64}, clustt2::PyArray{Int64}, minimum_excess_index::Int64)

    clustt1 = pyconvert(Array{Int64}, clustt1) 
    clustt2 = pyconvert(Array{Int64}, clustt2) 
    minimum_excess_index = pyconvert(Int64, minimum_excess_index)

    # Get lengths of cardinalities + switch around clusterings (such that an excess of workers is possible, not an excess of jobs) 
    lenn1, lenn2 = length(Set(clustt1)), length(Set(clustt2))

    # Flip them around if required (if the set of 2 is larger than the set of 1)
    if lenn2 > lenn1
        clust1 = clustt2 
        clust2 = clustt1 
    else
        clust1 = clustt1 
        clust2 = clustt2 
    end 

    # Contingency matrix 
    freqmat = freqtable(clust1, clust2)
    freqmat = convert(Array{Int64}, freqmat)
    max_cost = max(freqmat...) 
    cost = max_cost .- freqmat 
    #println("cost", cost)

    # Hungarian algorithm + Get all the workers/jobs rows/columns 
    assignments = solve_hungarian(cost) 
    workers = sort!(collect(Set(clust1)))
    jobs = sort!(collect(Set(clust2)))
    total_workers = length(workers) 
    worker_indices = 1:total_workers 

    # Find the useful (non-negative assignments) - row indices of relabel interest 
    useful_assignments_indices = findall(assignments.!=0)
    #println("useful assignments", useful_assignments_indices)

    # Get the workers corresponding to these useful indices (which will replace the assignments)
    row_indd = getindex(workers, getindex(worker_indices, useful_assignments_indices))
    #println("row indices", row_indd)

    # Get the assignments given to these indices- the jobs) 
    col_indd = getindex(jobs, getindex(assignments, useful_assignments_indices))
    #println("col indices", col_indd) 

    # Get the number of rows and the length of the clustering set 
    pairlength = length(row_indd) 
    clustlength = length(clust2)

    # Take back to_zero as necessary 
    if lenn2 > lenn1
        row_ind, col_ind = col_indd, row_indd
        clust1 = clustt1
        clust2 = clustt2 
    else
        row_ind, col_ind = row_indd, col_indd 
        clust1 = clustt1
        clust2 = clustt2 
    end 

    #println("row indices", row_ind)
    #println("col indices", col_ind) 
    #println(clust2)

    # Generate the new relabelled clustering
    newclust = zeros(Int64, clustlength)
    @simd for i in 1:pairlength 
        @simd for j in 1:clustlength
            @inbounds if clust2[j] == col_ind[i]
                newclust[j] = row_ind[i]
            end 
        end 
    end 

    # Identify which clustering has a higher cardinality 
    if maximum(clust2) > maximum(clust1)

        # Get uniques 
        uniques = sort!(collect(Set(clust2)))

        # Find uniques not remapped 
        not_remapped = setdiff(uniques, col_ind)
        #println("notmap", not_remapped)
        not_remapped_replacement = zeros(Int64, length(not_remapped))

        # Get new indices for not_remapped 
        @simd for i in 1:length(not_remapped)
            @inbounds not_remapped_replacement[i] = minimum_excess_index + i 
        end 

        # Relabel our newclust appropriately 
        @simd for i in 1:length(not_remapped)
            @simd for j in 1:clustlength
                @inbounds if clust2[j] == not_remapped[i] 
                    @inbounds newclust[j] = not_remapped_replacement[i] 
                end 
            end 
        end 


    end 
    
    return newclust 

end 
