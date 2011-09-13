require "kin_util"
require "generator"
ga_opts =  { verbose=1, population=100, generations=100, stop_fitenss=1e10,
      threads=5, seed_mul=3, sigfpe=true,
      converge=true, minpack = { maxfev = 1e6 }} 
filesave = "save.lua"
model = { 2, 2, 2 }
b     = #model-1
mp    = 2
mv    = 2
ma    = 1
-- Actual data for manipulation
P = {}
for i=1,mp do P[i] = {} for j=1,b do P[i][j] = rand_G() end end
if mv > 0 then
   V  = {}
   for i=1,mv do V[i] = {} for j=1,b do V[i][j] = merge( rand_plucker() ) end end
end
if ma > 0 then
   A  = {}
   for i=1,ma do A[i] = {} for j=1,b do A[i][j] = merge( rand_plucker() ) end end
end
-- Generate the model
out = io.open( "/tmp/model.lua", "w" )
generator_input( out, model, P, V, A )
out:close()
-- Load the model
dofile( "/tmp/model.lua" )
s   = gen_syn()
s:printClaim()
s:printJacobian( 1e-1 )


