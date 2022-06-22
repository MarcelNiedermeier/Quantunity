
#####################################
# Plotting script for random circuits
#####################################

using DelimitedFiles
using Plots
#using PlotlyJS
#using LaTeXStrings
gr()
#plotlyjs()
#pyplot()


# get data
datafile = readdlm("../Data/Random_circuits/random_circuit_N_12_bond_16_samples_100.csv", ',')

# general information ("header")
N = Int(datafile[1, 1])
maxdim = Int(datafile[1, 2])
bitstring = Int(datafile[1, 3])
length_depths = Int(datafile[1, 4])
N_samples = Int(datafile[1, 5])
depths = zeros(length_depths)
for i in 1:length_depths
    depths[i] = datafile[1, 5+i]
end

# data
PT_prob = datafile[2, :]
cum_distributions = zeros(length_depths, N_samples)
for i in 1:length_depths
    cum_distributions[i, :] = datafile[2+i, :]
end

# x for cumulative distribution
max_prob = 7
x_prob = LinRange(0, max_prob, N_samples)

# plotting
p1 = plot(x_prob, PT_prob, label="Porter-Thomas")
for i in 1:length(depths)
    plot!(x_prob, cum_distributions[i, :], label="N = $(N), χ = $(maxdim), D = $(depths[i])")
end
xlabel!("2^N  ρ")
ylabel!("P(p_U(x) <= ρ)")

savefig("../Plots/Random_circuits/random_circuit_N_$(N)_bond_$(maxdim)_samples_$(N_samples).png")
display(p1)
