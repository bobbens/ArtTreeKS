--[[
   Misc. nice useful Lua functions.
--]]


function merge(...)
    local t = {}
    for n = 1,select("#",...) do
        local arg = select(n,...)
        if type(arg)=="table" then
            for _,v in ipairs(arg) do
                t[#t+1] = v
            end
        else
            t[#t+1] = arg
        end
    end
    return t
end

function map(func, array)
   local new_array = {}
   for i,v in ipairs(array) do
      new_array[i] = func(v)
   end
   return new_array
end

function curry2(f)
   return function (x) return function (y) return f(x,y) end end
end


