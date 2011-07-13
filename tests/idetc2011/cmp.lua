

require "dir_path"

fname       = { "saved_syn", "gen_syn" }


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


function syn_cmp( file_a, file_b )
   syn_a = load_syn( file_a )
   syn_b = load_syn( file_b )
   return syn_a == syn_b
end


-- Code by David Kastrup
require "lfs"

function dirtree(dir)
  assert(dir and dir ~= "", "directory parameter is missing or empty")
  if string.sub(dir, -1) == "/" then
    dir=string.sub(dir, 1, -2)
  end

  local function yieldtree(dir)
    for entry in lfs.dir(dir) do
      if entry ~= "." and entry ~= ".." then
        entry=dir.."/"..entry
   local attr=lfs.attributes(entry)
   coroutine.yield(entry,attr)
   if attr.mode == "directory" then
     yieldtree(entry)
   end
      end
    end
  end

  return coroutine.wrap(function() yieldtree(dir) end)
end


list = {}
for filename, iter in dirtree( dir_path ) do
   list[ #list+1 ] = filename
end

eq_syn = {}
for i=1,#list-1 do
   for j=i+1,#list do
      local equal = syn_cmp( list[i], list[j] )
      if equal then
         eq_syn[ #eq_syn+1 ] = { list[i], list[j] }
      end
      print( equal, list[i], list[j] )
   end
end


if #eq_syn > 0 then
   print( "EQUAL:" )
   for k,v in ipairs(eq_syn) do
      print( tostring(k)..":   "..v[1].." - "..v[2] )
   end
end
print( string.format(
[[RESULTS:
   %d Analyzed
   %d Equal]], #list, #eq_syn ) )




