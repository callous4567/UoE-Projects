module fast_ising 

# The object that will be passed to functions to be iterated (has the matrix, size, everything.)
mutable struct twod_ising

    _mat::Array{Int64, 2}
    _M::Int64 
    _E::Float64
    _lx::Int64 
    _T::Float64 
    _dyn::String
    _binrate::Int64
    _max_sweeps::Int64
    _current_sweep::Int64
    _identifier::Int64 
    
end

# Periodic Boundary Conditions
function pbc(_i::Int64, _lx::Int64)

    if _i > (_lx)
        return Int64(1) 
    elseif _i < 1
        return Int64(_lx)
    else 
        return _i 
    end 

end 

# Magnetisation calculation
function fast_mag(_mat::Array{Int64,2})

    return sum(_mat)

end

# Energy calculation
function fast_energy(_mat::Array{Int64,2})

    _E = 0 
    _size = size(_mat)
    _lx,_ly = _size 

    for i in collect(1:_lx)

        for j in collect(1:_ly)

            _nnsum = _mat[pbc(i-1,_lx),j] +
                _mat[pbc(i+1,_lx),j] +
                _mat[i, pbc(j-1,_ly)] +
                _mat[i, pbc(j+1,_ly)] 
            _E += -1*_mat[i,j]*_nnsum 

        end 

    end 

    return Float64(_E/2)

end 

# Re-written. Use _ for local variables, no _ for global variables, etc) : Step the matrix forward by one
@fastmath function fast_glauber(_twod::twod_ising)

    # Define lx/etc 
    _lx = _twod._lx 

    # Generate random point 
    _i, _j = rand(1:_lx)[1], rand(1:_lx)[1]

    # Nearest neighbour sum 
    _nnsum = _twod._mat[pbc(_i-1,_lx),_j] + 
        _twod._mat[pbc(_i+1,_lx),_j] + 
        _twod._mat[_i, pbc(_j-1,_lx)] + 
        _twod._mat[_i, pbc(_j+1,_lx)]

    # Magnetisation cost of the change 
    _mag_cost = -2 * _twod._mat[_i,_j] 

    # Energy cost in terms of this 
    _energy_cost::Float64 = -1 * _mag_cost * _nnsum 

    # If/else of whether to do the change or not (accept if <= 0 else do probability check.)
    if _energy_cost <= 0

        _twod._mat[_i,_j] *= -1 
        _twod._M += _mag_cost 
        _twod._E += _energy_cost 

    else 

        _p = exp(-1 * _energy_cost / _twod._T)

        if rand(1)[1] <= _p

            _twod._mat[_i,_j] *= -1 
            _twod._M += _mag_cost 
            _twod._E += _energy_cost 

        end

    end 

    # Sweep complete 
    _twod._current_sweep += 1 

    return _twod 

end 

end 