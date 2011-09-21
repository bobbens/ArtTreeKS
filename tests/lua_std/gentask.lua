

require "generator"
require "kin_util"
require "util"


model = { 0, 3 }
b     = #model-1
mp    = 3
mv    = 2
ma    = 0


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
out = io.open( "model.lua", "w" )
generator_input( out, model, P, V, A )
out:close()

dofile( "model.lua" )
s = gen_syn()

