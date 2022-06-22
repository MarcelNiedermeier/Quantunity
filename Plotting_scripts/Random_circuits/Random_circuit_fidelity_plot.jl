

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

calculate_ex = false
random = true

# get data
datafile = readdlm("Data/Random_circuits/random_circuit_fid_D_100_bond_32_samples_20.csv", ',')
if calculate_ex
    datafile_ex = readdlm("Data/Random_circuits/random_circuit_fid_ex_D_80_bond_16_samples_20.csv", ',')
end
if random
    datafile_rand = readdlm("Data/Random_circuits/random_circuit_rand_MPS_fid_D_100_bond_32_samples_20.csv", ',')
end

# general information ("header")
N_samples = Int(datafile[1, 1])
maxdim = Int(datafile[1, 2])
bitstring = Int(datafile[1, 3])
length_Ns = Int(datafile[1, 4])
depth = Int(datafile[1, 5])
Ns = zeros(length_Ns)
for i in 1:length_Ns
    Ns[i] = Int(datafile[1, 5+i])
end

# data
twofidelities = zeros(length_Ns, depth÷2+1)
Nfidelities = zeros(length_Ns, depth÷2+1)
Avtwofidelities = zeros(length_Ns, depth÷2+1)
twofidelities_std = zeros(length_Ns, depth÷2+1)
Nfidelities_std = zeros(length_Ns, depth÷2+1)
Avtwofidelities_std = zeros(length_Ns, depth÷2+1)

twofidelities_rand = zeros(length_Ns, depth÷2+1)
Nfidelities_rand = zeros(length_Ns, depth÷2+1)
Avtwofidelities_rand = zeros(length_Ns, depth÷2+1)

twofidelities_ex = zeros(length_Ns, depth÷2+1)
Nfidelities_ex = zeros(length_Ns, depth÷2+1)
Avtwofidelities_ex = zeros(length_Ns, depth÷2+1)

for i in 1:length_Ns
    twofidelities[i, :] = datafile[1+i, :]
    Nfidelities[i, :] = datafile[1+length_Ns+i, :]
    Avtwofidelities[i, :] = datafile[1+2*length_Ns+i, :]
    twofidelities_std[i, :] = datafile[1+3*length_Ns+i, :]
    Nfidelities_std[i, :] = datafile[1+4*length_Ns+i, :]
    Avtwofidelities_std[i, :] = datafile[1+5*length_Ns+i, :]
end

if random
    for i in 1:length_Ns
        twofidelities_rand[i, :] = datafile_rand[1+i, :]
        Nfidelities_rand[i, :] = datafile_rand[1+length_Ns+i, :]
        Avtwofidelities_rand[i, :] = datafile_rand[1+2*length_Ns+i, :]
    end
end

if calculate_ex
    for i in 1:length_Ns
        twofidelities_ex[i, :] = datafile_ex[i, :]
        Nfidelities_ex[i, :] = datafile_ex[length_Ns+i, :]
        Avtwofidelities_ex[i, :] = datafile_ex[2*length_Ns+i, :]
    end
end


##########
# plotting
##########

if calculate_ex
    for i in 1:length_Ns
        println(twofidelities[i, :] - twofidelities_ex[i, :])
    end
end


depths = collect(1:depth÷2+1)

labels = reshape([Symbol("N = $(Int(Ns[i]))") for i in 1:length(Ns)], 1, length(Ns))
labels2 = reshape([Symbol("N = $(Int(Ns[i])), geom. av.") for i in 1:length(Ns)], 1, length(Ns))
labels3 = reshape([Symbol("N = $(Int(Ns[i])), rand. in. state") for i in 1:length(Ns)], 1, length(Ns))

p1 = plot(depths, transpose(twofidelities), lab=map(string, labels))
plot!(depths, transpose(Avtwofidelities), lab=map(string, labels2))
xlabel!("D")
ylabel!("2-qubit fidelity")
display(p1)
savefig("Plots/Random_circuits/random_circuit_fid_D_$(depth)_bond_$(maxdim)_samples_$(N_samples).png")

if calculate_ex
    p2 = plot(depths, transpose(twofidelities_ex), lab=map(string, labels))
    plot!(depths, transpose(Avtwofidelities_ex), lab=map(string, labels2))
    xlabel!("D")
    ylabel!("2-qubit fidelity")
    display(p2)
    #savefig("Plots/Random_circuits/random_circuit_fid_D_$(depth)_bond_$(maxdim)_samples_$(N_samples).png")
end

if random
    p3 = plot(depths, transpose(twofidelities), lab=map(string, labels))
    plot!(depths, transpose(twofidelities_rand), lab=map(string, labels3))
    xlabel!("D")
    ylabel!("2-qubit fidelity")
    display(p3)
end

#p3 = plot(depths, abs.(transpose(twofidelities)-transpose(twofidelities_ex)), lab=map(string, labels), yaxis=:log)
#plot!(depths, transpose(Avtwofidelities_ex), lab=map(string, labels2))
#xlabel!("D")
#ylabel!("2-qubit fidelity")
#display(p3)



#p2 = plot(depths, transpose(Avtwofidelities), lab=map(string, labels))
#xlabel!("D")
#ylabel!("av. 2-qubit fidelity")
#display(p2)
#savefig("Plots/Random_circuits/random_circuit_av_fid_D_$(depth)_bond_$(maxdim)_samples_$(N_samples).png")
