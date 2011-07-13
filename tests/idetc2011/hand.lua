

require "util"
require "kin_util"
require "luadq"


function table_chop( t, ts, te )
   local nt = {}
   for i=ts,te do
      nt[ #nt+1 ] = t[i]
   end
   return nt
end

require "hand_data"
common   = table_chop( joint,  1,  3 )
index    = table_chop( joint,  4,  7 )
middle   = table_chop( joint,  8, 11 )
third    = table_chop( joint, 12, 16 )
fourth   = table_chop( joint, 17, 21 )
thumb    = table_chop( joint, 22, 26 )


hand = {}
hand.axes = { common, index, middle, third, fourth, thumb }
hand.fingers = table_chop( hand.axes, 2, #hand.axes )

hand.G = {}
for i=1,5 do
   hand.G[i] = luadq.homo( R[i], T[i] )
   hand.G[i] = luadq.inv( hand.G[i] )
end

hand.branches = {}
for i=1,#hand.fingers do
   hand.branches[i]  = merge( common, hand.fingers[i] )
end

hand.tcp = function ( angles )
   local P = {}
   for i=1,#hand.fingers do
      local thetas = merge( angles[1], angles[1+i] )
      P[i] = kin_fk( thetas, hand.branches[i], hand.G[i] )
   end
   return P
end




