--------------------------------------------------------------------------------
--- Examples:
--- time th ./test.lua -s 1 -e 4000
--------------------------------------------------------------------------------
local torch = require('torch')
local moses = require('moses')
local bio = require('bio')
local xlua = require('xlua')

local cmd = torch.CmdLine()
cmd:option('-i', 'allU.fasta', 'fasta file input')
cmd:option('-d', 'multithreaded-rnagrep', 'output directory')
cmd:option('-o', 'temp-testout.csv', 'output file name')
cmd:option('-p', 'crumple', 'Name of program to run, default is crumple')
cmd:option('-w', '30', 'window size')
cmd:option('-z', '1', 'which sequence to parse')
cmd:option('-s', '1', 'start position for slider')
cmd:option('-e', '100', 'end position for slider')
local opt = cmd:parse(arg or {})

local method = ("echo $winseq | $method"):gsub('$method', opt.p)

local seqs = bio.read(opt.i)
for i, j in ipairs(seqs) do
   seqs[i].win = bio.win(seqs[i].na, tonumber(opt.w))
end

if tonumber(opt.e) > #seqs[1].win then
   opt.e = #seqs[1].win
end

local cleanup = function(t,pos,last)
   local temp = {}
   if (pos == last) then;
      for j = 1, #t do
	 temp[j] = t[j]
      end
   else
      for j = 1, 25 do
	 temp[j] = t[(25+j)]
      end
   end
   return temp
end

-- for z = 1, #seqs do
local z = tonumber(opt.z)
-- print(seqs[z].name)
local mfile = torch.MemoryFile()
local g_hp = {}
local s_hp = {}
local e_hp = {}
local w_hp = {}
local inputs = tonumber(opt.s)
local inpute = tonumber(opt.e)
for i = inputs, inpute do
   -- print("DOING  " .. seqs[z].name .. "  " .. i .. "/" .. opt.e)
   xlua.progress(i, inpute) 
   if (i%50) == 0  then;
      s_hp = cleanup(s_hp,i,inpute)
      e_hp = cleanup(e_hp,i,inpute)
      w_hp = cleanup(w_hp,i,inpute)
      g_hp = cleanup(g_hp,i,inpute)
   end
   
   local cmd = method:gsub('$winseq', seqs[z].win[i])
   local f = assert(io.popen(cmd, 'r'))
   for line in f:lines() do
      local s, e, g  = string.find(line, "(%b())")
      
      if (g ~= nil) and
	 (string.match(g, "%(%.%.%.%.?%.?%.?%.?%.?%.?%)") ~= nil) and
	 (string.match(g, "%(%.%.%.%.%.%.%.%.%.%.+%)") == nil) and
	 (string.match(g, "%(.*%).*%(") == nil) and
	 (string.match(g, "%).*%(.*%)") == nil) and
	 (string.match(g, "%(%.%.%.%.%.+%(") == nil) and 
	 (string.match(g, "%)%.%.%.%.%.+%)") == nil) and
	 (string.match(g, "%(+%.+%(+%.+%(+") == nil) and
	 (string.match(g, "%)+%.+%)+%.+%)+") == nil) and
	 (string.match(g, "%(%(%(%(+") ~= nil) and
	 (string.match(g, "%)%)%)%)+") ~= nil)
      then;
	 
	 local store = true
	 for j = 1, #g_hp do
	    
	    if ((((i+s-1)==s_hp[j]) and ((i+e-1)==e_hp[j])) and ((g==g_hp[j])))
	    then;
	       store = false
	       break
	    end
	    
	 end

	 if (store == true) then;
	    s_hp[#s_hp+1] = (i+s-1)
	    e_hp[#e_hp+1] = (i+e-1)
	    g_hp[#g_hp+1] = g
	    w_hp[#w_hp+1] = string.sub(seqs[z].win[i],s,e)
	    mfile:writeString(s_hp[#s_hp] .. "," .. e_hp[#s_hp] .. "," .. w_hp[#s_hp] .. "," .. g_hp[#s_hp] .. "\n")
	 end
	 
      end
   end
   f:close()
   
end

mfile:seek(1)
local df = mfile:readString("*a")
mfile:close()
-- df = table.concat(df)
local fout = assert(io.open("./" .. opt.d .. "/" .. seqs[z].name .. "_" .. opt.s  .. "_" .. opt.e .. "_" .. opt.o .. ".csv", "w"))
fout:write(df)
fout:close()
