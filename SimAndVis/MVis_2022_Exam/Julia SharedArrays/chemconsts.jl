@everywhere const dx = Float64(1)
@everywhere const dt = Float64(0.1)
@everywhere const lx = Int64(256)
@everywhere const ly = Int64(256)
@everywhere const lz = Int64(256)
@everywhere const shape = (lx,ly,lz)
@everywhere const area = lx*ly 
@everywhere const length = lx*ly*lz 
@everywhere const D = Float64(0.5)
@everywhere const q = Float64(2.5)
@everywhere const p = Float64(0.5)
@everywhere const max_sweeps = Int64(10000)
@everywhere const binrate = Int64(10)
@everywhere const c1 = D*dt/(dx^2)
@everywhere const c2 = q*dt 
@everywhere const c3 = p*dt 
@everywhere const indices = (1,2,3,4)