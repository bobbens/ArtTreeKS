
require "dir_path"

out      = dir_path.."/exp%03d.lua"
data_path = dir_path.."/data.csv"
echo     = "echo "
echo     = ""

-- Detect starting i
i  = 0
repeat
   i  = i + 1
   local s = string.format( out, i )
   if fp ~= nil then
      fp:close()
   end
   fp = io.open( s )
until fp == nil
print( "Starting at "..tostring(i) )

ga_opts =  { verbose=1, population=100, generations=100, stop_fitenss=1e10,
   threads=5, seed_mul=3, sigfpe=true,
   converge=true, minpack = { maxfev = 1e6 }} 

function process_file( filename, filesave, datout )
   dofile( filename )
   local s = gen_syn()
   local sts = s:stats()
   print( string.format( "n=%d, m=%d, ni=%d, mi=%d, L=%d, r=%d, b=%d", sts.n, sts.m, sts.ni, sts.mi, sts.L, sts.r, sts.b ) )
   local res = s:solver_ga( ga_opts )
   local e = res.elapsed
   print( string.format( "fitness = %.3e in %d:%02d:%02d", res.fit_best, e/3600, (e/60)%60, e%60) )
   print( "Saving to "..filesave )
   s:save( filesave )

   -- Open output
   local dout
   local exists = false
   dout     = io.open( datout, "r" )
   if dout ~= nil then
      exists = true
      dout:close()
   end
   dout     = io.open( datout, "a" )
   if not exists then
      dout:write( "i,best,generations,elapsed\n" )
   end
   dout:write( string.format( "%d,%e,%d,%d", i, res.fit_best, res.generations, res.elapsed ) )
   dout:close()
end

-- To infinity and beyond
while true do
   local sout = string.format( out, i )
   process_file( "test_run.lua", sout, data_path )
   i = i + 1
end


