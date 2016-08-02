--------------------------------------------------------------------------------
--- Examples:
--- time th /files/projects/RNAgrep/v1.0/compare.lua -input " .csv, .csv, .csv, .csv, .csv, .csv" -nseq "wt,pthree,twothree,wtwo,lansing,cox" -output "natural.csv"
--- time th /files/projects/RNAgrep/v1.0/compare.lua -input " .csv, .csv, .csv" -nseq "wt,max,sd" -output "synthetic.csv"
--- time th /files/projects/RNAgrep/v1.0/compare.lua -input '/files/projects/RNAgrep/v1.0/filtered/wt_1_6603_v1.csv,/files/projects/RNAgrep/v1.0/filtered/23127_1_6603_v1.csv,/files/projects/RNAgrep/v1.0/filtered/cox_1_6603_v1.csv,/files/projects/RNAgrep/v1.0/filtered/lansing_1_6603_v1.csv,/files/projects/RNAgrep/v1.0/filtered/max_1_6603_v1.csv,/files/projects/RNAgrep/v1.0/filtered/P3L37_1_6603_v1.csv,/files/projects/RNAgrep/v1.0/filtered/sd_1_6603_v1.csv,/files/projects/RNAgrep/v1.0/filtered/w2_1_6603_v1.csv' -nseq "wt,twothree,cox,lansing,max,pthree,sd,wtwo" -output "compared-all.csv"
--------------------------------------------------------------------------------
local torch = require('torch')
local csv = require('ftcsv')
local xlua = require('xlua')
local re = require('re')
--- ++++++ Arguements parse
local cmd = torch.CmdLine()
cmd:option('-input', 'input.csv', 'comma separated file names of input files')
cmd:option('-nseq', 'seq,names', 'comma seperated names for each input file')
cmd:option('-output', 'output.csv', 'output file')
local opt = cmd:parse(arg or {})


--- ++++++ Load csv 
local input = csv.parse(opt.input .. "\r\n", ",", {loadFromString = true, headers = false})
local t = {}
local csv_options = {headers = false, rename = {"start", "stop", "seq", "hp"}}
for k, v in ipairs(input[1]) do
   t[k]  = csv.parse(v, ",", csv_options)
end

--- ++++++ Load text file containing each seqs name, in correct order
local seqnames = csv.parse(opt.nseq .. "\r\n", ",", {loadFromString = true, headers = false})
pat = re.compile[[ 
listname <- {| ({[a-z][a-z]*} %s* ','*)* |} 
]]
local seqnames = pat:match(opt.nseq)
print(seqnames)
-- local f = assert(io.open(opt.seqnames, 'r'))
-- for line in f:lines() do
--    seqnames[#seqnames+1] = line
-- end

--- ++++++ dumping to memory
local mfile = torch.MemoryFile()

--- ++++++ Make comparison with first seq as the reference
local matches = {}
local count   = tonumber(#t[1])
for i = 1, count do
   xlua.progress(i, count)
   matches[i] = {}
   matches[i][seqnames[1]] = {}
   matches[i][seqnames[1]].start    = t[1][i].start 
   matches[i][seqnames[1]].stop     = t[1][i].stop
   matches[i][seqnames[1]].seq      = t[1][i].seq
   matches[i][seqnames[1]].hp       = t[1][i].hp   
   for k = 2, #seqnames do
      matches[i][seqnames[k]] = {}
      matches[i][seqnames[k]].start    = "NA"
      matches[i][seqnames[k]].stop     = "NA"
      matches[i][seqnames[k]].seq      = "NA"
      matches[i][seqnames[k]].hp       = "NA"   
      local kcount = tonumber(#t[k])
      for j = 1, kcount do
	 if (
	    (t[1][i].start       ==    t[k][j].start) and
	       (t[1][i].stop     ==    t[k][j].stop ) and 
	       (t[1][i].hp       ==    t[k][j].hp   )
	 )
	 then
	    matches[i][seqnames[k]].start    = t[k][j].start 
	    matches[i][seqnames[k]].stop     = t[k][j].stop
	    matches[i][seqnames[k]].seq      = t[k][j].seq
	    matches[i][seqnames[k]].hp       = t[k][j].hp
	 end
      end
   end
end

local c = ""
for i = 1, #seqnames do
   if i == #seqnames then 
      mfile:writeString(c .. seqnames[i] .. "-start" .."," .. seqnames[i] .. "-stop" .."," .. seqnames[i] .. "-seq" .."," .. seqnames[i] .. "-hp" .. "\n")
      c = ""
   else
      c = c .. seqnames[i] .. "-start" .."," .. seqnames[i] .. "-stop" .."," .. seqnames[i] .. "-seq" .."," .. seqnames[i] .. "-hp" .. ", "
   end
end

c = ""
for i = 1, #matches do
   for j = 1, #seqnames do
      if j == #seqnames then 
	 mfile:writeString(c .. matches[i][seqnames[j]].start .."," .. matches[i][seqnames[j]].stop .."," .. matches[i][seqnames[j]].seq .."," .. matches[i][seqnames[j]].hp .. "\n")
	 c = ""
      else
	 c = c .. matches[i][seqnames[j]].start .."," .. matches[i][seqnames[j]].stop .."," .. matches[i][seqnames[j]].seq .."," .. matches[i][seqnames[j]].hp .. ", "
      end
   end
end

mfile:seek(1)
local df = mfile:readString("*a")
mfile:close()
local file = assert(io.open(opt.output, "w"))
file:write(df)
file:close()

