-- Checks the output of two PairHMM implementations for correctness
--
-- User: carneiro
-- Date: 3/13/13

if table.maxn(arg) < 2 then
    print("You must provide the two files to compare")
    print("lua checkLikelihoods.lua <expected.out> <actual.out>")
    os.exit(1);
end

local function getValues(e, a)
    return e:read("*n"), a:read("*n")
end

local expected = io.open(arg[1])
local actual   = io.open(arg[2])
local sensitivity = arg[3] or 0.001

local e, a = getValues(expected, actual)
local fail = 0
local total = 0
while (e ~= nil and a ~= nil) do
    if math.abs(a - e) > sensitivity then
        print(total .. ":", a, e, "["..math.abs(a-e).."]")
        fail = fail + 1
    end
    total = total + 1
    e, a = getValues(expected, actual)
end

if fail == 0 then
    print("all " .. total .. " tests passed.")
else
    print(fail .. " out of " .. total .. " tests failed.")
end