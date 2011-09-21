
require "luadq"
require "prng"


function kin_fk( thetas, axes, G )
   local P = luadq.point( { 0, 0, 0 } )
   local Q = luadq.point( { 0, 0, 0 } )
   for i=1,#axes do
      Q = Q*luadq.rotation_plucker( thetas[i], axes[i][1], axes[i][2] )
   end
   Q = Q*G
   P = P:f4g( Q )
   return P
end

function dot( u, v )
   return u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
end

function cross( u, v )
   return { u[2]*v[3] - u[3]*v[2],
           -u[1]*v[3] + u[3]*v[1],
            u[1]*v[2] - u[2]*v[1] }
end

function frames_needed( b, r )
   return math.ceil((4*r) / (6*b - r)) + 1
end

function rand_G ()
   local R, T
   local s, s0 = rand_plucker()
   local a = prng.num() * math.pi * 2
   local d = prng.num() * 10
   R = luadq.rotation( a, s, s0 )
   T = luadq.translation( d, s )
   return R * T
end

function rand_vec_n ()
   local s = { prng.num(), prng.num(), prng.num() }
   local m = math.sqrt( dot(s, s) )
   for i=1,3 do s[i] = s[i] / m end
   return s
end

function rand_plucker ()
   local s, c, s0
   s = rand_vec_n()
   c = { prng.num(), prng.num(), prng.num() }
   for i=1,3 do c[i] = c[i]*10 end
   s0 = cross( c, s )
   return s, s0
end

function rand_angles( r )
   local angles = {}
   for i=1,r do
      angles[i] = prng.num() * math.pi * 2
   end 
   return angles
end


