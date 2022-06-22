

#####################################################
# Plotting script for QFT fidelity, compared to exact
#####################################################

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
datafile = readdlm("../Data/Quantum_Fourier_Transform/QFT_MPS_exact_fid_maxN_16_maxbond_80_samples_4.csv", ',')





D = Int(datafile[1, 1])
N_sample = Int(datafile[1, 2])
len_Ns = Int(datafile[1, 3])
Ns = datafile[2, 1:1+len_Ns-1]
maxdims = map(Int, datafile[3, :])

av_fidelities = zeros(len_Ns, length(maxdims))
for i in 1:len_Ns
    av_fidelities[i, :] = datafile[3+i, :]
end

#println([Symbol("χ0 = $(Int(randombonds[i]))") for i in 1:length(randombonds)])

# plotting
labels = reshape([Symbol("N = $(Int(Ns[i]))") for i in 1:length(Ns)], 1, length(Ns))
#p1 = plot(maxdims, av_fidelities[1, :], lab=map(string, labels))
#for i in 2:len_randombonds
#    plot!(maxdims, av_fidelities[i, :], lab=map(string, labels))
#end
p1 = plot(maxdims, transpose(av_fidelities), lab=map(string, labels), legend=:bottomright)
xlabel!("bond dim")
ylabel!("fidelity")
# QFT_MPS_fid_N_40_maxbond_64_samples_4.
savefig("../Plots/Quantum_Fourier_Transform/QFT_MPS_fid_exact_maxN_$(Ns[end])_maxbond_$(maxdims[end])_samples_$(N_sample).png")
display(p1)
