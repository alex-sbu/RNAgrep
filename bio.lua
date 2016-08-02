local torch = require('torch')
local re = require('re')
local ftcsv = require('ftcsv')
local bio = {}


bio = {

   fasta = (
      function(f)
	 local file = io.open(f, "r")
	 local seqname = ""
	 local seqna = {}
	 local open = false
	 local seq = {}
	 local i = 0
	 for line in file:lines() do	    
	    if string.match(line, ">") == nil and open == true then
	       table.insert(seqna, line)
	    end	    
	    if #seqna >= 1 then
	       seq[i] = {
		  name        = seqname,
		  na          = table.concat(seqna)
	       }
	    end
	    if string.match(line, ">") ~= nil and open == true then
	       seqname = ""
	       seqna = {}
	       open = false
	    end	    
	    if string.match(line, ">") ~= nil and open == false then
	       seqname = string.gsub(line, ">(.+)", "%1")
	       open = true
	       i = i+1
	    end
	 end
	 return seq
      end
	  ),

   flatten = (
      --- ++++++ source: http://svn.wildfiregames.com/public/ps/trunk/build/premake/premake4/src/base/table.lua
      function(arr)
	 local result = {}	 
	 local function toflat(arr)
	    for _, v in ipairs(arr) do
	       if type(v) == "table" then
		  toflat(v)
	       else
		  table.insert(result, v)
	       end
	    end
	 end	 
	 toflat(arr)
	 return result
      end
	     ),
   
   select = (
      --- ++++++ source: http://stackoverflow.com/questions/31304601/torch-lua-how-to-select-a-subset-of-an-array-or-tensor
      function(t, first, last)
	 local sub = {}
	 for i=first,last do
	    sub[#sub + 1] = t[i]
	 end
	 return sub
      end      
	    ),
   
   win = (
      function(seq, winsize)
	 local seql = string.len(seq)-winsize
	 local seqw = {}
	 local seqr = torch.range(1, seql)
	 seqr:apply(function(wstart)
	       local wend = wstart+winsize-1
	       seqw[wstart] = string.sub(seq, wstart, wend)
		    end
	 )
	 return seqw
      end
	 ),

   transpose = (
      function(t)
	 local tb = {}
	 for row = 1, #t[1] do
	    tb[row] = {}
	    for col = 1, #t do
	       tb[row][col] = t[col][row]
	    end
	 end
	 return tb
      end
	       ),

   csv = (
      function(inp)
	 local hotdatain = nil
	 if(inp.data == nil) then	                       -- READ IN
	    if (inp.headers == nil) then 
	       inp.headers = false
	    end
	    if (string.match(inp.file, ",") ~= nil) then      -- MANY FILES
	       hotdatain = {}
	       local files = {}
	       local files = string.split(inp.file, ",")
	       for k, v in ipairs(files) do
		  hotdatain[k]  = bio.csvin(v, inp.headers)
	       end
	    else		                              -- SINGLE FILE
	       hotdatain = bio.csvin(inp.file, inp.headers)
	    end
	 else                                                 -- WRITE OUT
	    if (torch.type(inp.headers) == "string") then
	       bio.csvout(inp.data, inp.file, inp.headers)    -- ASSIGN HEADERS
	    else
	       bio.csvout(inp.data, inp.file)                 -- NO HEADERS
	    end
	 end
	 if (hotdatain ~= nil) then
	    return hotdatain
	 end
      end
	 ),

   csvin = (
      function(file, headers)
	 local csvopt = { headers = headers }
	 local datain = ftcsv.parse(file, ",", csvopt)
	 return datain
      end
	   ),

   csvout = (
      function(data, saveas, headers)	 
	 if string.match(saveas, ".+%.csv") ~= nil then
	    saveas = saveas .. ".csv"
	 end
	 if string.match(torch.type(data), "torch.+") ~= nil then
	    if (headers ~= nil) then
	       bio.tensor2csv(data, saveas, headers)          -- TENSOR TO CSV 
	    else
	       bio.tensor2csv(data, saveas)
	    end
	 elseif torch.type(data) == "table" then 
	    local file = assert(io.open(saveas, "w"))   
	    if (headers ~= nil) then
	       file:write(headers .. "\n")
	    end
	    file:write(ftcsv.encode(data, ","))              -- LUA TABLE TO CSV
	    file:close()
	 end
      end      
	    ),
   
   tensor2csv = (
      function(m, f, headers)
	 local ne = 0
	 local drillm = m
	 local c = 0
	 while string.match(torch.type(data), "torch.+") ~= nil do
	    print(torch.type(drillm))
	    c = c + 1
	    if torch.type(drillm) ~= "number" then
	    drillm = drillm[1]
	    end
	    if string.match(torch.type(data), "torch.+") == nil then
	       ne = c
	    end
	 end
	 if ne == 3 then
	    for x = 1, m:size()[1] do
	       local t = {}
	       for i = 1, m:size()[2] do
		  for j = 1, m:size()[3] do
		     if j ~= m:size()[3] then
			t[(#t+1)] = m[x][i][j] .. ","
		     else
			t[(#t+1)] = m[x][i][j] .. "\n"
		     end
		  end
	       end
	       file = assert(io.open(x .. "-" .. f, "w"))
	       if (headers ~= nil) then
		  file:write(headers .. "\n")
	       end
	       file:write(table.concat(t, ""))
	       file:close()
	    end
	 end
	 if ne == 2 then
	    local t = {}
	    for i = 1, m:size()[1] do
	       for j = 1, m:size()[2] do
		  if j ~= m:size()[2] then
		     t[(#t+1)] = m[i][j] .. ","
		  else
		     t[(#t+1)] = m[i][j] .. "\n"
		  end
	       end
	    end 
	    file = assert(io.open(f, "w"))
	    if (headers ~= nil) then
	       file:write(headers .. "\n")
	    end
	    file:write(table.concat(t, ""))
	    file:close()
	 end
	 if ne == 1 then
	    local t = {}
	    for i = 1, m:size()[1] do
	       t[(#t+1)] = m[i] .. "\n"
	    end
	    file = assert(io.open(f, "w"))
	    if (headers ~= nil) then
	       file:write(headers .. "\n")
	    end
	    file:write(table.concat(t, ""))
	    file:close()
	 end
      end
)
}

return bio
