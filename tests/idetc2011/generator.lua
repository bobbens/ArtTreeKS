
require "luadq"
require "prng"
require "util"
require "kin_util"


bounds   = false
funccall = false

-- Define curry stuff
function add(x,y) return x+y end
function sub(x,y) return y-x end
cadd = curry2( add )
csub = curry2( sub )

-- Handle parameters
seed = 9
if seed < 0 then
   math.randomseed( tonumber(tostring(os.time()):reverse():sub(1,6)) )
   seed = math.random( 1e6 )
end
prng.init( seed )


function dof_calc_r( dof )
   if type(dof) == "table" then
      local r = 0
      for i=1,#dof do
         r = r + dof_calc_r( dof[i] )
      end
      return r
   elseif type(dof) == "number" then
      return dof
   end
   return 0
end


function dof_calc_b( dof )
   if type(dof) == "table" then
      local b = 0
      for i=1,#dof do
         b = b + dof_calc_b( dof )
      end
      return b
   elseif type(dof) == "number" then
      return 0
   elseif dof == "t" then
      return 1
   end
   return 0;
end


function create_chain( name, out, axes, angles, info, m )
   out:write( "local "..name.." = kin_object.new( \"chain\" )\n" )
   for j=1,#axes do
      -- Create joint
      out:write( "local j = kin_joint.new( \"revolute\" )\n" )

      -- Write joint axis
      local s  = axes[j][1]
      local s0 = axes[j][2]
      out:write( string.format( "j:setPlucker( { %.18e, %.18e, %.18e }, { %.18e, %.18e, %.18e } ) -- dot=%.3e\n",
            s[1], s[2], s[3], s0[1], s0[2], s0[3], dot(s, s0)  ) )

      -- Write axis bounds
      if info.bounds then
         local val = 10
         out:write( string.format( "j:setPluckerBounds( { %.18e, %.18e, %.18e }, { %.18e, %.18e, %.18e },\n"
                                .. "                    { %.18e, %.18e, %.18e }, { %.18e, %.18e, %.18e } )\n",
               -1, -1, -1, 1, 1, 1, s0[1]-val, s0[2]-val, s0[3]-val, s0[1]+val, s0[2]+val, s0[3]+val ) )
      end

      -- Write joint angles
      out:write( "j:setPositions( { " )
      for i=2,m-1 do
         out:write( tostring( angles[j][i] - angles[j][1] ) .. ", " )
      end
      out:write( tostring( angles[j][m] - angles[j][1] ) .. " } )\n" )

      -- Write angle bounds
      if info.bounds then
         local lb = ""
         local ub = ""
         for i=2,m do
            local c
            if i < m then c = ", " else c = " " end
            lb = lb .. string.format( "%.18e"..c, 0 )
            ub = ub .. string.format( "%.18e"..c, 2*math.pi )
         end
         out:write( "j:setPositionBounds( { "..lb.." },\n"
                 .. "                     { "..ub.." } )\n" )
      end

      -- Attach
      out:write( name..":attach( j )\n" )
   end
end


function create_tcp( name, out, P, info, m )
   out:write( "local "..name.." = kin_object.new( \"tcp\" )\n" )
   out:write( name..":setFK( {\n" )
   for j=1,m do
      local R,d = P[j]:extract()
      out:write( "   -- Frame " .. tostring(j) .. "\n" )
      out:write( "   {\n" )
      for i=1,3 do
         out:write( string.format( " { %.18e, %.18e, %.18e, %.18e },\n", R[i][1], R[i][2], R[i][3], d[i] ) )
         out:write( "    " )
      end
      out:write( " { 0., 0., 0., 1. } " )
      if j==m then out:write( "}\n" ) else out:write( "},\n" ) end
   end
   out:write( "   } )\n" )
end


function create_splitter( name, out, info )
   out:write( "local "..name.." = kin_object.new( \"splitter\" )\n" )
end


function create_n_axes( n, info )
   local t = {}
   for i=1,n do
      t[ #t+1 ] = { { 1, 0, 0 }, { 0, 0, 0 } }
   end
   return t
end


function vsum( v )
   local s = 0
   for _,i in ipairs(v) do
      s = s + i
   end
   return s
end

function table_chop( t, ts, te )
   local nt = {}
   for i=ts,te do
      nt[ #nt+1 ] = t[i]
   end
   return nt
end

function generator_from( out, n, axes, in_angles, G, P )
   local info        = {}
   info.bounds       = true
   local angles      = {}
   local b           = #n-1
   local r           = vsum( n )
   local m           = frames_needed( b, r )
   local name        = string.format("%d-{",n[1])
   for i=1,b do
      if i ~= 1 then
         name = name..","
      end
      name = name..string.format("%d",n[1+i])
   end
   name = name.."}"

   for i=1,#n do -- Axis
      angles[i]   = {}
      for j=1,n[i] do -- Chain
         angles[i][j] = {}
         for k=1,m do -- Pose
            angles[i][j][k] = in_angles[k][i][j]
         end
      end
   end

   out:write( "require \"synthesis\"\n" )
   out:write( "function gen_syn ()\n" )
   out:write( "local s = syn.new( "..tostring(m).." )\n" )
   create_splitter( "spl",  out,  info )
   create_chain(    "ko",  out,  axes[1], angles[1], info, m )
   for i=1,b do
      local ko  = string.format( "ko%d",  i )
      local tcp = string.format( "tcp%d", i )
      create_chain( ko, out, axes[1+i], angles[1+i], info, m )
      local Pv  = {}
      for j=1,m do
         Pv[j] = P[j][i]
      end
      create_tcp(  tcp, out, Pv, info, m )
      out:write( ko..":attach( "..tcp.." )\n" )
      out:write( "spl:attach( "..ko.." )\n" )
   end
   out:write( "ko:attach( spl )\n" )
   out:write( "s:addObject( ko )\n" )
   out:write( "s:finalize()\n" )
   out:write( "return s\n" )
   out:write( "end\n" )
end


function generator_input( out, n, P )
   local info        = {}
   info.bounds       = true
   local axes        = {}
   local angles      = {}
   local b           = #n-1
   local r           = vsum( n )
   local m           = frames_needed( b, r )
   local name        = string.format("%d-{",n[1])
   for i=1,b do
      if i ~= 1 then
         name = name..","
      end
      name = name..string.format("%d",n[1+i])
   end
   name = name.."}"
   for i=1,#n do -- Axis
      axes[i]     = create_n_axes( n[i] )
      angles[i]   = {}
      for j=1,n[i] do -- Chain
         angles[i][j] = {}
         for k=1,m do -- Pose
            angles[i][j][k] = 0
         end
      end
   end

   out:write( "require \"synthesis\"\n" )
   out:write( "function gen_syn ()\n" )
   out:write( "local s = syn.new( "..tostring(m).." )\n" )
   create_splitter( "spl",  out,  info )
   create_chain(    "ko",  out,  axes[1], angles[1], info, m )
   for i=1,b do
      local ko  = string.format( "ko%d",  i )
      local tcp = string.format( "tcp%d", i )
      create_chain( ko, out, axes[1+i], angles[1+i], info, m )
      local Pv  = {}
      for j=1,m do
         Pv[j] = P[j][i]
      end
      create_tcp(  tcp, out, Pv, info, m )
      out:write( ko..":attach( "..tcp.." )\n" )
      out:write( "spl:attach( "..ko.." )\n" )
   end
   out:write( "ko:attach( spl )\n" )
   out:write( "s:addObject( ko )\n" )
   out:write( "s:finalize()\n" )
   out:write( "return s\n" )
   out:write( "end\n" )
end

function tbl_str( t )
   local s = "{ "
   for i=1,#t do
      if i ~= 1 then
         s = s..", "
      end
      s = s..tostring(t[i])
   end
   return s.." }"
end


