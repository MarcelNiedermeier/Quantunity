
using ITensors

# check single site gates as MPO
xtest = OpSum()
xtest += "Sz*Sz", 5
#xtest += "ProjUp", 5
#xtest += "ProjDn", 5

# get sites
sites = siteinds("Qubit",N)

# turn into MPO
gateX = MPO(xtest, sites)

# has bond dimensions of 3?
println("checking X gate: ", gateX)
