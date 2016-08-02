bk = require('bk-tree')
torch = require('torch')
xlua = require('xlua')
bio = require('bio')
moses = require('moses')

local cmd = torch.CmdLine()
cmd:option('-f', 'lev.csv', 'input file')
cmd:option('-d', 'compared-all', 'temporary output directory for data processing')
local opt = cmd:parse(arg or {})


dstr   = bio.csv({file = "./".. opt.d .."/structures.csv"})
dstart = bio.csv({file = "./".. opt.d .."/start.csv"})
dend   = bio.csv({file = "./".. opt.d .."/end.csv"})
table.remove(dstr)
table.remove(dstart)
table.remove(dend)
dstr   = bio.flatten(dstr)
dstart = bio.flatten(dstart)
dend   = bio.flatten(dend)
duni = moses.unique(dstr)
du = bio.select(duni, 1, #duni)
dr = bio.select(dstr, 1, #dstr)
ds = bio.select(dstart, 1, #dstart)
de = bio.select(dend, 1, #dend)

print("Making Levenshtein distance matrix")
tm = torch.Tensor(tonumber(#dr), tonumber(#du)):zero()
for k, v in ipairs(du) do
   xlua.progress(k, #du)
   for q, z in ipairs(dr) do
      tm[{q, k}] = bk.levenshtein_dist(v, z)
   end
end

t = {}
for i = 1, tm:size()[1] do
   for j = 1, tm:size()[2] do
      if j == tm:size()[2] then
	 t[(#t+1)] = tm[{i,j}]  .. "\n"
      else
	 t[(#t+1)] = tm[{i,j}]  .. ","
      end
   end
end
file = assert(io.open("/files/projects/RNAgrep/v1.0/src/".. opt.d .."/" .. opt.f, "w"))
file:write(table.concat(t, ""))
file:close()

