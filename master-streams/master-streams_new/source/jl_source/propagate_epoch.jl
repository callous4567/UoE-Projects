using PythonCall
mastorad = 1/206264806.2471

"""
This aims to be a Julia port of the PyGaia Epoch Propagation code. 
Any and all credits and licenses for this code lie with the original author- I merely ported it.
Use at your own risk. 
"""

# Number of seconds in Julian year
julian_year_seconds = 365.25 * 86400.0
# Astronomical Unit in meter, IAU constant and defining length
au_in_meter = 149597870700.0
# AU expressed in km*yr/s
au_km_year_per_sec = au_in_meter / (julian_year_seconds * 1000.0)

@fastmath function normal_triad(
    phi::Float64, theta::Float64
    )
    p = Vector{Float64}(undef, 3)
    q = Vector{Float64}(undef, 3)
    r = Vector{Float64}(undef, 3)
    @inbounds p[1],p[2],p[3] = -1*sin(phi),cos(phi),0
    @inbounds q[1],q[2],q[3] = -1*sin(theta)*cos(phi),-1*sin(theta)*sin(phi),cos(theta)
    @inbounds r[1],r[2],r[3] = cos(theta)*cos(phi),cos(theta)*sin(phi),sin(theta)
    return p,q,r 
end 

@fastmath function cartesian_to_spherical(x::Float64, y::Float64, z::Float64)

    rCylSq = x^2 + y^2 
    r = sqrt(rCylSq + z^2)
    phi = atan(y, x) # replace with atan if precision is required. 
    return r, phi, atan(z, sqrt(rCylSq))

end 

@fastmath function _propagate_astrometry(
    phi::Float64, theta::Float64, parallax::Float64, muphistar::Float64, mutheta::Float64, vrad::Float64, t0::Int64, t1::Int64
    )

    t = t1 - t0
    p0, q0, r0 = normal_triad(phi, theta) # vectors 

    pmra0 = muphistar*mastorad
    pmdec0 = mutheta*mastorad
    pmr0 = vrad*parallax / au_km_year_per_sec*mastorad 

    pmtot0sqr = (muphistar^2 + mutheta^2)*mastorad^2

    pmvec0 = pmra0.*p0 + pmdec0.*q0 # vector 

    f = (1 + 2*pmr0*t + (pmtot0sqr + pmr0^2)*t^2) ^ (-0.5)
    u = (r0*(1 + pmr0*t) .+ pmvec0.*t).*f # vector 

    _, phi1, theta1 = cartesian_to_spherical(u[1], u[2], u[3])
    parallax1 = parallax*f
    pmr1 = (pmr0 + (pmtot0sqr + pmr0^2) * t)*f^2
    pmvec1 = (pmvec0.*(1 + pmr0*t) - r0.*(pmr0^2 * t)).* f^3 # vector 
    p1, q1, r1 = normal_triad(phi1, theta1) # vector 
    muphistar1 = sum(p1.*pmvec1./mastorad) # vector 
    mutheta1 = sum(q1.*pmvec1./mastorad) # vector 
    murad1 = pmr1/mastorad 

    # Convert phi1 to degrees and make sure its [0,360] 
    if phi1 < 0
        phi1 += 2*pi 
    end 
    phi1, theta1 = rad2deg(phi1), rad2deg(theta)

    return phi1, theta1, parallax1, muphistar1, mutheta1, murad1

end 

@fastmath function propagate_astrometry(
    phi::PyArray{Float64}, theta::PyArray{Float64}, parallax::PyArray{Float64}, muphistar::PyArray{Float64}, mutheta::PyArray{Float64}, vrad::PyArray{Float64}, t0::Int64, t1::Int64
    )

    """
    Same as the PyGaia documentation. Propagate the astrometric parameters of a source from the reference epoch t0 to
    the new epoch t1. Note coordinates are ra-dec/lat-long style, not theta-phi style. 

    Parameters
    ----------
    phi : float
        Longitude at reference epoch (radians).
    theta : float
        Latitude at reference epoch (radians). 
    parallax : float
        Parallax at the reference epoch (mas).
    muphistar : float
        Proper motion in longitude (including np.cos(latitude) term) at reference
        epoch (mas/yr).
    mutheta : float
        Proper motion in latitude at reference epoch (mas/yr).
    vrad : float
        Radial velocity at reference epoch (km/s).
    t0 : float
        Reference epoch (Julian years).
    t1 : float
        New epoch (Julian years).

    Returns
    -------
    phi1, theta1, parallax1, muphistar1, mutheta1, murad1 : float or array
        Astrometric parameters, including the "radial proper motion" (NOT the radial
        velocity), at the new epoch.
    """

    phi = pyconvert(Vector{Float64}, phi)
    theta = pyconvert(Vector{Float64}, theta)
    parallax = pyconvert(Vector{Float64}, parallax)
    muphistar = pyconvert(Vector{Float64}, muphistar)
    mutheta = pyconvert(Vector{Float64}, mutheta)
    vrad = pyconvert(Vector{Float64}, vrad)

    phi1, theta1, parallax1, muphistar1, mutheta1, murad1 = Vector{Float64}(undef, size(phi)[1]),Vector{Float64}(undef, size(phi)[1]),Vector{Float64}(undef, size(phi)[1]), 
                                                            Vector{Float64}(undef, size(phi)[1]),Vector{Float64}(undef, size(phi)[1]),Vector{Float64}(undef, size(phi)[1])
    @simd for i in 1:size(phi)[1]
        @inbounds phi1[i], theta1[i], parallax1[i], muphistar1[i], mutheta1[i], murad1[i] = 
        _propagate_astrometry(phi[i], theta[i], parallax[i], muphistar[i], mutheta[i], vrad[i], t0, t1)
    end 

    return phi1, theta1, parallax1, muphistar1, mutheta1, murad1 

end 

