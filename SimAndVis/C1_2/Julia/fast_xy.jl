module fast_xy 

using Tullio
using Distributions 
using LinearAlgebra

mutable struct xy_ising

    _mat::Array{Float64, 3}
    _M::Vector{Float64}
    _E::Float64
    _lx::Int64 
    _T::Float64 
    _binrate::Int64 
    _max_sweeps::Int64 
    _current_sweep::Int64
    _identifer::Int64 

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

# Mag calculation
@fastmath function fast_mag(_mat::Array{Float64, 3})

    _mag = Vector{Float64}(undef, 2)
    @tullio _mag[1] = _mat[i,j,1]
    @tullio _mag[2] = _mat[i,j,2]

    return _mag 
            
end 

# Energy calculation
@fastmath function fast_energy(_mat::Array{Float64, 3})

    _E::Float64 = 0
    _size = size(_mat)

    for _j in 1:_size[2]

        for _i in 1:_size[1] 

            _E += dot(_mat[_i,_j,:], 
            _mat[pbc(_i-1,_size[1]),_j,:] +
            _mat[pbc(_i+1,_size[1]),_j,:] +
            _mat[_i, pbc(_j-1,_size[2]),:] +
            _mat[_i, pbc(_j+1,_size[2]),:])

        end 

    end 

    return _E 

end 

# Evaluate atan2 fully-manually 
@fastmath function atan2(y::Float64, x::Float64)
    if x > 0
        return atan(y/x)
    elseif x < 0 
        if y >= 0
            return atan(y/x) + pi
        else 
            return atan(y/x) - pi 
        end 
    else 
        if y > 0 
            return pi/2 
        elseif y < 0
            return -1*pi/2 
        else 
            return 0
        end 
    end 
end 

# Approximate method (njuffa!)
# https://math.stackexchange.com/questions/1098487/atan2-faster-approximation
@fastmath function atan2_approx(y::Float64, x::Float64)
    a = min(abs(x), abs(y))/max(abs(x), abs(y))
    s = a^2 
    r = ((-0.0464964749*s + 0.15931422)*s - 0.327622764)*s*a + a
    if abs(y) > abs(x)
        r = 1.57079637 - r
    end 
    if x < 0
        r = 3.14159274 - r
    end 
    if y < 0
        r = -r
    end 
    return r 
end 


# Angles calculation (0,2pi)
@fastmath function _angle(_mat::Array{Float64, 3})


    _size = size(_mat) 

    __angle = Array{Float64}(undef, _size[1:2])

    for _j in 1:_size[2]

        for _i in 1:_size[1]

            __angle[_i,_j] = atan2_approx(_mat[_i,_j,2], _mat[_i,_j,1])

            if __angle[_i,_j] < 0
                __angle[_i,_j] += 2*pi

            end 

        end 

    end  

    return __angle 

# With no approx

end 
    
# Generate a random normalized normal in 2D 
function autoij() 
    _i, _j = rand(Normal(0,1)), rand(Normal(0,1))
    _mag = sqrt(_i^2 + _j^2)
    _i /= _mag 
    _j /= _mag 
    _vect = Vector{Float64}(undef, 2)
    _vect[1], _vect[2] = _i,_j 
    return _vect 
end 

# Autoij but without returning vectors 
function autoij_novect() 
    _i, _j = rand(Normal(0,1)), rand(Normal(0,1))
    _mag = sqrt(_i^2 + _j^2)
    _i /= _mag 
    _j /= _mag 
    return _i,_j 
end 

# Autoangle but without approximations 
@fastmath function noapprox_autoangle()
    _i, _j = rand(Normal(0,1)), rand(Normal(0,1))
    _mag = sqrt(_i^2 + _j^2)
    _i /= _mag 
    _j /= _mag 
    uuu = atan(_j,_i)
    if uuu < 0
        uuu += 2*pi 
    end 
    return uuu
end 


# Re-written. Use _ for local variables, no _ for global variables, etc) : Step the matrix forward by one
@fastmath function fast_glauber(_mat::Array{Float64,3}, _lx::Int64, _M::Vector{Float64}, _E::Float64, _T::Float64)

    # Generate random point 
    _i, _j = rand(1:_lx), rand(1:_lx)

    # Nearest neighbour sum 
    _nnsum = _mat[pbc(_i-1,_lx),_j,:] + 
        _mat[pbc(_i+1,_lx),_j,:] + 
        _mat[_i, pbc(_j-1,_lx),:] + 
        _mat[_i, pbc(_j+1,_lx),:]

    # Select a new angle for the spin (i,j) (unit vector)
    _new_ij = autoij()

    # Get magnetisation cost 
    _mag_cost = _new_ij - _mat[_i,_j,:]

    # Energy cost in terms of this 
    _energy_cost = -1 * dot(_mag_cost, _nnsum)

    # If/else of whether to do the change or not (accept if <= 0 else do probability check.)
    if _energy_cost <= 0

        _mat[_i,_j, :] = _new_ij
        _M += _mag_cost 
        _E += _energy_cost 

    else 

        _p = exp(-1 * _energy_cost / _T)

        if rand() <= _p

            _mat[_i,_j,:] = _new_ij 
            _M += _mag_cost 
            _E += _energy_cost 

        end

    end 

    return _mat, _M, _E 

end 

# Work with angles directly. Note that using an approx tan here won't work due to the nature of MCMC
@fastmath function fast_angleglauber(_mat::Array{Float64,2}, _lx::Int64, _T::Float64)

    # Generate random point 
    _i, _j = rand(1:_lx), rand(1:_lx)

    # Select a new angle for the spin (i,j) (unit vector)
    _new_ij = noapprox_autoangle()

    # Evaluate the change in angle 
    _dij = _new_ij - _mat[_i,_j]

    # Energy cost in terms of this 
    _energy_cost = -2*(sin(_mat[pbc(_i-1,_lx),_j] - _mat[_i,_j] - _dij/2) +
        sin(_mat[pbc(_i+1,_lx),_j] - _mat[_i,_j] - _dij/2) +
        sin(_mat[_i,pbc(_j-1,_lx)] - _mat[_i,_j] - _dij/2) +
        sin(_mat[_i,pbc(_j+1,_lx)] - _mat[_i,_j] - _dij/2))*sin(_dij/2)
    

    # If/else of whether to do the change or not (accept if <= 0 else do probability check.)
    if _energy_cost <= 0

        _mat[_i,_j] = _new_ij 

    else 

         _p = exp(-1 * _energy_cost / _T)

        if rand() <= _p

            _mat[_i,_j] = _new_ij 

        end

    end 

    return _mat

end 


# convert oned to twod indices for square grids 
function oned_twod(_lx::Int64, __i::Int64)
    
    # Get the _i index 
    _i = floor(Int64, (__i-1)/_lx)

    # Next the _j index
    _j = floor(Int64, __i - _i*_lx)

    return _i+1, _j 
end 

# convert twod to oned indices for square grids. NO BOUNDS CHECK!
function twod_oned(_lx::Int64, _i::Int64, _j::Int64)

    return (_i - 1)*_lx + _j 

end 


# Work with angles directly. Note that using an approx tan here won't work due to the nature of MCMC
@fastmath function fast_onedangleglauber(_vec::Vector{Float64}, _lx::Int64, _T::Float64)

    # Generate random point (contiguous random) 
    _i, _j = rand(1:_lx), rand(1:_lx)

    # Select a new angle for the spin (i,j) (unit vector)
    _new_ij = noapprox_autoangle()

    # Evaluate the change in angle 
    _dij = _new_ij - _vec[twod_oned(_lx,_i,_j)]

    # Energy cost in terms of this 
    _energy_cost = -2*(sin(_vec[twod_oned(_lx, pbc(_i-1,_lx),_j)] - _vec[twod_oned(_lx,_i,_j)] - _dij/2) +
        sin(_vec[twod_oned(_lx,pbc(_i+1,_lx),_j)] - _vec[twod_oned(_lx,_i,_j)] - _dij/2) +
        sin(_vec[twod_oned(_lx,_i,pbc(_j-1,_lx))] - _vec[twod_oned(_lx,_i,_j)] - _dij/2) +
        sin(_vec[twod_oned(_lx,_i,pbc(_j+1,_lx))] - _vec[twod_oned(_lx,_i,_j)] - _dij/2))*sin(_dij/2)
    

    # If/else of whether to do the change or not (accept if <= 0 else do probability check.)
    if _energy_cost <= 0

        _vec[twod_oned(_lx,_i,_j)] = _new_ij 

    else 

         _p = exp(-1 * _energy_cost / _T)

        if rand() <= _p

            _vec[twod_oned(_lx,_i,_j)] = _new_ij 

        end

    end 

    return _vec

end 













end 

