

fname = { "saved_syn", "gen_syn" }


function load_syn( filename )
   if filename == nil then
      error( "Filename is nil!" )
   end

   for _,v in ipairs(fname) do
      _G[v] = nil
   end

   dofile( filename )

   for _,v in ipairs(fname) do
      if _G[v] ~= nil then
         return _G[v]()
      end
   end

   error( "Unable to load syn from "..filename )
end


s  = load_syn( arg[1] )
if arg[2] ~= nil then
   s2 = load_syn( arg[2] )
end
print( arg[2] )
--s:print()
if s2 ~= nil then
   s:visualize( s2 )
else
   s:visualize()
end

