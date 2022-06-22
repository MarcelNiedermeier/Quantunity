

####################################
# Plotting script for QFT (fidelity)
####################################

using DelimitedFiles
using Plots
#using PlotlyJS
#using LaTeXStrings
gr()
#plotlyjs()
#pyplot()

#calculate_ex = false
#random = true

# get data
datafile = readdlm("../Data/Quantum_Fourier_Transform/QFT_MPS_fid_N_8_maxbond_28_samples_4.csv", ',')

N = Int(datafile[1, 1])
N_sample = Int(datafile[1, 2])
len_randombonds = Int(datafile[1, 3])
randombonds = datafile[2, 1:1+len_randombonds-1]
maxdims = map(Int, datafile[3, :])

av_fidelities = zeros(len_randombonds, length(maxdims))
for i in 1:len_randombonds
    av_fidelities[i, :] = datafile[3+i, :]
end

println([Symbol("χ0 = $(Int(randombonds[i]))") for i in 1:length(randombonds)])

# plotting
labels = reshape([Symbol("χ0 = $(Int(randombonds[i]))") for i in 1:length(randombonds)], 1, length(randombonds))
#p1 = plot(maxdims, av_fidelities[1, :], lab=map(string, labels))
#for i in 2:len_randombonds
#    plot!(maxdims, av_fidelities[i, :], lab=map(string, labels))
#end
p1 = plot(maxdims, transpose(av_fidelities), lab=map(string, labels), legend=:bottomright)
xlabel!("bond dim")
ylabel!("fidelity")
# QFT_MPS_fid_N_40_maxbond_64_samples_4.
savefig("../Plots/Quantum_Fourier_Transform/QFT_MPS_fid_N_$(N)_maxbond_$(maxdims[end])_samples_$(N_sample).png")
display(p1)
