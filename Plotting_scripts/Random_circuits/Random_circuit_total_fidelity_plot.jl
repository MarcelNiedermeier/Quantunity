


################################################
# Plotting script for random circuits (fidelity)
################################################

using DelimitedFiles
using Plots
#using PlotlyJS
#using LaTeXStrings
gr()
#plotlyjs()
#pyplot()


# get data
datafile = readdlm("Data/Random_circuits/random_circuit_fid_D_100_N_60_samples_3.csv", ',')

# general information ("header")
N_samples = Int(datafile[1, 1])
N = Int(datafile[1, 2])
bitstring = Int(datafile[1, 3])
length_bond_dims = Int(datafile[1, 4])
depth = Int(datafile[1, 5])
bond_dims = zeros(length_bond_dims)
for i in 1:length_bond_dims
    bond_dims[i] = Int(datafile[1, 5+i])
end

# data
twofidelities = zeros(length_bond_dims, depth÷2+1)
Nfidelities = zeros(length_bond_dims, depth÷2+1)
Avtwofidelities = zeros(length_bond_dims, depth÷2+1)
twofidelities_std = zeros(length_bond_dims, depth÷2+1)
Nfidelities_std = zeros(length_bond_dims, depth÷2+1)
Avtwofidelities_std = zeros(length_bond_dims, depth÷2+1)

for i in 1:length_bond_dims
    twofidelities[i, :] = datafile[1+i, :]
    Nfidelities[i, :] = datafile[1+length_bond_dims+i, :]
    Avtwofidelities[i, :] = datafile[1+2*length_bond_dims+i, :]
    twofidelities_std[i, :] = datafile[1+3*length_bond_dims+i, :]
    Nfidelities_std[i, :] = datafile[1+4*length_bond_dims+i, :]
    Avtwofidelities_std[i, :] = datafile[1+5*length_bond_dims+i, :]
end

Avtwofidelities_final = zeros(length_bond_dims)
for i in 1:length_bond_dims
    Avtwofidelities_final[i] = Avtwofidelities[i, end]
end

Avtwoerror = zeros(length_bond_dims)
for i in 1:length_bond_dims
    Avtwoerror[i] = 1 - Avtwofidelities_final[i]
end


##########
# plotting
##########

depths = collect(1:depth÷2+1)

labels = reshape([Symbol("χ = $(Int(bond_dims[i]))") for i in 1:length(bond_dims)], 1, length(bond_dims))
labels2 = reshape([Symbol("χ = $(Int(bond_dims[i])), geom. av.") for i in 1:length(bond_dims)], 1, length(bond_dims))
labels3 = reshape([Symbol("N = $N")], 1, 1)

p1 = plot(depths, transpose(twofidelities), lab=map(string, labels))
plot!(depths, transpose(Avtwofidelities), lab=map(string, labels2))
xlabel!("D")
ylabel!("2-qubit fidelity")
display(p1)
savefig("Plots/Random_circuits/random_circuit_fid_D_$(depth)_N_$(N)_samples_$(N_samples).png")

p2 = plot(depths, transpose(Nfidelities), lab=map(string, labels))#, yaxis=:log)
xlabel!("D")
ylabel!("N-qubit fidelity")
display(p2)
savefig("Plots/Random_circuits/random_circuit_N_fid_D_$(depth)_N_$(N)_samples_$(N_samples).png")

p3 = plot(bond_dims, Avtwofidelities_final, lab=map(string, labels))#, xaxis=:log)
xlabel!("χ")
ylabel!("av. 2-qubit fidelity")
display(p3)
savefig("Plots/Random_circuits/random_circuit_av_fid_D_$(depth)_N_$(N)_samples_$(N_samples).png")

p4 = plot(bond_dims, Avtwoerror, lab=map(string, labels3))#, xaxis=:log)
xlabel!("χ")
ylabel!("av. 2-qubit error")
display(p4)
savefig("Plots/Random_circuits/random_circuit_av_error_D_$(depth)_N_$(N)_samples_$(N_samples).png")
