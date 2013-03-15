--
-- Created by IntelliJ IDEA.
-- User: carneiro
-- Date: 3/13/13
-- Time: 11:44 PM
-- To change this template use File | Settings | File Templates.
--

if table.maxn(arg) < 3 then
    print("You must provide the test data and the evaluation results")
    print("lua checkHaplotypes.lua <input_data> <expected_likelihoods> <actual_likelihoods>")
    os.exit(-1)
end

-- command line arguments
local testData = io.open(arg[1])
local expected = io.open(arg[2])
local actual = io.open(arg[3])

-- find maximum likelihood haplotype given field (actual or expected)
local function findMax (t, field) 
	local max = nil
	local maxHap = nil
	for h, v in pairs(t) do
		if max == nil or v[field] > max then 
			maxHap = h
		   	max = v[field]
		end
	end
	return maxHap
end

-- facility function that returns both the actual and expected max likelihood haplotypes
local function checkMostLikelyHaplotypes(t)
    local a, e = findMax(t, "actual"), findMax(t, "expected")
	local matches, mismatches = 0, 0
	if (a == e) then 
		matches = 1
	else 
		print("MISMATCH")
		print("actual: ".. a)
		print("expected: " .. e)
		print(t[a].actual, t[e].actual)
		print(t[a].expected, t[e].expected)
		mismatches = 1
	end
	return matches, mismatches
end


local function initLEA(pl, pe, pa) 
	local l = pl or testData:read()
	local e = pe or expected:read("*n")
	local a = pa or actual:read("*n")
	return l, e, a 
end

-- find all haplotypes while reading the first block
-- todo: special case a block with only 1 read and 1 haplotype
local function findAllHaplotypes(hapHash, pl, pa, pe)
	local l, e, a = initLEA(pl, pe, pa)
	local hap, read = l:match("(%w+)%s(%w+)")
	while hapHash[hap] == nil do
		hapHash[hap] = {hap = hap, actual = a, expected = e}
		l = testData:read()
		hap, read = l:match("(%w+)%s(%w+)")
		a, e = actual:read("*n"), expected:read("*n")
	end
	return l, e, a
end

local function addLikelihoods(hapHash, pl, pe, pa)
	local l, e, a = initLEA(pl, pe, pa)
	local hap = l:match("(%w+)%s(%w+)")
	while hapHash[hap] ~= nil do
		hapHash[hap].actual = hapHash[hap].actual + a
		hapHash[hap].expected = hapHash[hap].expected + e
		l = testData:read()
		hap, read = l:match("(%w+)%s(%w+)")
		a, e = actual:read("*n"), expected:read("*n")
	end
	return l, e, a
end

local function dump (x)
	for i, t in pairs(x) do
		print(i..":")
		for k, v in pairs(t) do
			print("\t", k, v)
		end
	end
end

local hapHash = {}
local l, e, a 
local i = 1
local matches, mismatches = 0, 0
repeat 
	l, e, a = findAllHaplotypes(hapHash, l, e, a)
	if l == nil then break end -- special case for truncated inputs
	l, e, a = addLikelihoods(hapHash, l, e, a)
	local m, mm = checkMostLikelyHaplotypes(hapHash)
	matches = matches + m
	mismatches = mismatches + mm
	hapHash = {}
	i = i + 1 
	if mm > 0 then print(matches, mismatches) end
until l == nil 

print("Total matches: " .. matches)
print("Total mismatches: " .. mismatches)


