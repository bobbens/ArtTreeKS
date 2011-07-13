

require "util"
require "hand"
require "generator"
require "prng"


angles_range = { -- Common (Wrist)
                 { 0, math.pi/4 }, -- Z
                 { 0, math.pi/4 }, -- Y
                 { 0, math.pi/4 }, -- X
                 -- Index
                 { 0, math.pi/4 }, -- Y  CMC
                 { 0, -math.pi/4 }, -- Z  MCP
                 { 0, math.pi/2.2 }, -- Y  PIP
                 { 0, math.pi/4.5 }, -- Y  DIP
                 -- Middle
                 { 0, math.pi/4 }, -- Y  CMC
                 { 0, -math.pi/8 }, -- Z  MCP
                 { 0, math.pi/2 }, -- Y  PIP
                 { 0, math.pi/4 }, -- Y  DIP
                 -- Third
                 { 0, math.pi/14 }, -- Y  CMC 1
                 { 0, math.pi/4 }, -- Y  CMC 2
                 { 0, math.pi/8 }, -- Z  MCP
                 { 0, math.pi/2 }, -- Y  PIP
                 { 0, math.pi/4 }, -- Y  DIP
                 -- Fourth
                 { 0, math.pi/11 }, -- Y  CMC 1
                 { 0, math.pi/4 }, -- Y  CMC 2
                 { 0, math.pi/4 }, -- Z  MCP
                 { 0, math.pi/2 }, -- Y  PIP
                 { 0, math.pi/4 }, -- Y  DIP
                 -- Thumb
                 { 0, math.pi/8 }, -- CMC 1
                 { 0, math.pi/6 }, -- CMC 2
                 { 0, -math.pi/3 },--math.pi/4 }, -- MCP 1
                 { 0, -math.pi/5 },---math.pi/4 }, -- MCP 2
                 { 0, math.pi/2}--math.pi/2 }  -- IP
                 }

-- We only really want 9 frames but we generate 33
mf = 33
m  = 9 -- overkill

P      = {}
angles = {}
for f=1,mf do -- All frames
   angles[f] = {}
   local n   = #hand.axes
   local acc = 1
   for i=1,n do -- All axes in the frame
      angles[f][i] = {}
      local nn     = #hand.axes[i]
      for j=1,nn do -- All joints in the axis

         -- Calculate basis
         local a = angles_range[acc]
         local r = a[2]-a[1]
         local val = a[1] + (f-1)/(m-1) * a[2]
         acc = acc + 1

         -- Add noise
         local nse = 0
         if f ~= 1 then
            nse = r/7 * prng.num()
         end

         -- Form
         angles[f][i][j] = val + nse
      end
   end

   P[f] = hand.tcp( angles[f] )
end

-- Here we can choose to generate for thumb or any other system of fingers
use_thumb = true
if use_thumb then
   haxes = { common, index, middle, thumb }
   hand.G[3]      = hand.G[5]
   for f=1,mf do
      P[f][3]        = P[f][5]
      angles[f][4]   = angles[f][6]
   end
else
   haxes = { common, index, middle, third }
end
n     = { 3, 3, 3, 3 }
--n     = { 3, 4, 4, 5 }
nh    = {} for i=1,#haxes do nh[ #nh+1 ] = #haxes[i] end
--[[
haxes = { common, index, middle, third, fourth, thumb }
nh    = {} for i=1,#haxes do nh[ #nh+1 ] = #haxes[i] end
n     = nh
--]]

-- Output files
out   = io.open( "test_run.lua",      "w" )
out2  = io.open( "test_run_base.lua", "w" )

-- Data generation
generator_from(   out2, nh, haxes, angles, hand.G, P )
generator_input(  out, n, P )

-- Close files
out:close()
out2:close()





