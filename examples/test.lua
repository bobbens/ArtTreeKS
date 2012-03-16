

require "util"
require "generator"
require "prng"


-- Options
tmp_name = "out_tmp.lua"
filesave = "out_result.lua"
ga_opts =  { verbose=1, population=100, generations=100, stop_fitenss=1e10,
      threads=5, seed_mul=3, sigfpe=true,
      converge=true, minpack = { maxfev = 1e6 }} 


function rand_vel ()
   local t = {}
   for i=1,6 do
      t[i] = math.pi*(prng.num()*2-1)
   end
   return t
end

-- Generate robot
n  = {3,2,2,2}
b  = #n-1
mp = 3
mv = 2
ma = 0
P  = {}
for i=1,mp do P[i] = {} for j=1,b do P[i][j] = rand_G() end end
if mv > 0 then
   V  = {}
   for i=1,mv do V[i] = {} for j=1,b do V[i][j] = rand_vel() end end
end
if ma > 0 then
   A  = {}
   for i=1,ma do A[i] = {} for j=1,b do A[i][j] = rand_vel() end end
end

-- Create output
out = io.open( tmp_name, "w" )
generator_input( out, n, P, V, A )
out:close()

-- Process
print( "---LOADING TEST---" )
dofile( tmp_name )
local s   = gen_syn()
--s:print()
s:solver_minpack()
--s:print()
local sts = s:stats()
print( string.format( "n=%d, m=%d, ni=%d, mi=%d, L=%d, r=%d, b=%d", sts.n, sts.m, sts.ni, sts.mi, sts.L, sts.r, sts.b ) )
local res = s:solver_ga( ga_opts )
local e   = res.elapsed
print( string.format( "fitness = %.3e in %d:%02d:%02d", res.fit_best, e/3600, (e/60)%60, e%60) )
s:print()
print( "Saving to "..filesave )
s:save( filesave )
--[[
--]]

-- Test save
if res ~= nil and res.fit_best > 1e10 then
   print( "---TESTING SAVE---" )
   dofile( filesave )
   local s   = saved_syn()
   local sts = s:stats()
   print( string.format( "n=%d, m=%d, ni=%d, mi=%d, L=%d, r=%d, b=%d", sts.n, sts.m, sts.ni, sts.mi, sts.L, sts.r, sts.b ) )
   s:visualize()
end


