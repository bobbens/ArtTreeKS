

filename = arg[1] or "2-2_3.lua"
filesave = "result.lua"


dofile( filename )
ga_opts =  { verbose=1, population=100, generations=100, stop_fitenss=1e10,
   threads=5, seed_mul=3, sigfpe=true,
   converge=true, minpack = { maxfev = 1e6 }} 
s = gen_syn()
--s:print()
--s:solver_minpack()
sts = s:stats()
print( string.format( "n=%d, m=%d, ni=%d, mi=%d, L=%d, r=%d, b=%d", sts.n, sts.m, sts.ni, sts.mi, sts.L, sts.r, sts.b ) )
res = s:solver_ga( ga_opts )
e = res.elapsed
print( string.format( "fitness = %.3e in %d:%02d:%02d", res.fit_best, e/3600, (e/60)%60, e%60) )
print( "Saving to "..filesave )
s:save( filesave )
--s:print()

