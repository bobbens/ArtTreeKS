

filename = arg[1] or error("File not specified")
filesave = "result.lua"


function opt_ga( s )
   local ga_opts = { verbose=1, population=100, generations=100, stop_fitenss=1e10,
      threads=5, seed_mul=3, sigfpe=true,
      converge=true, minpack = { maxfev = 1e6 }} 
   local res = s:solver_ga( ga_opts )
   local e = res.elapsed
   print( string.format( "fitness = %.3e in %d:%02d:%02d", res.fit_best, e/3600, (e/60)%60, e%60) )
   return true
end


function opt_cmaes( s )
   local cmaes_opts = {
         converge=false,
         lambda=1e4
         }
   local res = s:solver_cmaes( cmaes_opts )
   local e = res.elapsed
   print( string.format( "fitness = %.3e in %d:%02d:%02d", res.fit_best, e/3600, (e/60)%60, e%60) )
   return true
end


local opt_try = {
   --opt_ga,
   opt_cmaes
}

dofile( filename )
s = gen_syn()
--s:print()
--s:solver_minpack()
sts = s:stats()
print( string.format( "n=%d, m=%d, ni=%d, mi=%d, L=%d, r=%d, b=%d", sts.n, sts.m, sts.ni, sts.mi, sts.L, sts.r, sts.b ) )
for _,solver in ipairs(opt_try) do
   s = gen_syn()
   solver( s )
end
--print( "Saving to "..filesave )
--s:save( filesave )
--s:print()

