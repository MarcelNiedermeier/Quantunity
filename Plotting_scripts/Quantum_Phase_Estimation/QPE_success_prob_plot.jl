
##############################################
## Plotting script for QPE success probability
##############################################


using DelimitedFiles
using Plots
#using PlotlyJS
#using LaTeXStrings
gr()
#plotlyjs()
#pyplot()


# get data
datafile = readdlm("../Data/Quantum_Phase_Estimation/QPE_MPS_success_prob_n_prec_10_maxbond_40_samples_4.csv", ',')
datafile2 = readdlm("../Data/Quantum_Phase_Estimation/QPE_MPS_success_prob_2_n_prec_10_maxbond_40_samples_4.csv", ',')



# first data file
n_prec = Int(datafile[1, 1])
N_sample = Int(datafile[1, 2])
len_n_probs = Int(datafile[1, 3])
n_probs = datafile[2, 1:1+len_n_probs-1]
maxdims = map(Int, datafile[3, :])

av_success_probs = zeros(len_n_probs, length(maxdims))
for i in 1:len_n_probs
    av_success_probs[i, :] = datafile[3+i, :]
end


# second data file
n_prec2 = Int(datafile2[1, 1])
N_sample2 = Int(datafile2[1, 2])
len_maxdims = Int(datafile2[1, 3])
maxdims2 = datafile2[2, 1:1+len_maxdims-1]
n_probs2 = map(Int, datafile2[3, :])

av_success_probs2 = zeros(len_maxdims, length(n_probs2))
for i in 1:len_maxdims
    av_success_probs2[i, :] = datafile2[3+i, :]
end


""" Success probability for phase estimation with t-n = nprob
qubits. """
function success_prob(nprob)
    return 1 - 1/(2*exp(nprob) - 4)
end


#println([Symbol("χ0 = $(Int(randombonds[i]))") for i in 1:length(randombonds)])

# plotting
labels = reshape([Symbol("n_probs = $(Int(n_probs[i]))") for i in 1:length(n_probs)], 1, length(n_probs))
labels2 = reshape([Symbol("maxdim = $(Int(maxdims2[i]))") for i in 1:length(maxdims2)], 1, length(maxdims2))

#p1 = plot(maxdims, av_fidelities[1, :], lab=map(string, labels))
#for i in 2:len_randombonds
#    plot!(maxdims, av_fidelities[i, :], lab=map(string, labels))
#end

p1 = plot(maxdims, transpose(av_success_probs), lab=map(string, labels), legend=:bottomright)
xlabel!("bond dim")
ylabel!("success prob")
# QPE_MPS_success_prob_n_prec_6_maxbond_8_samples_4
savefig("../Plots/Quantum_Phase_Estimation/QPE_MPS_success_prob_n_prec_$(n_prec)_maxbond_$(maxdims[end])_samples_$(N_sample).png")
display(p1)

p2 = plot(n_probs2, transpose(av_success_probs2), lab=map(string, labels2), legend=:bottomright)
success_probs = [success_prob(i) for i in n_probs2]
plot!(n_probs2, success_probs, lab="analytical")
xlabel!("n prob")
ylabel!("success prob")
# QPE_MPS_success_prob_2_n_prec_6_maxbond_8_samples_4
savefig("../Plots/Quantum_Phase_Estimation/QPE_MPS_success_prob_2_n_prec_$(n_prec)_maxbond_$(maxdims[end])_samples_$(N_sample).png")
display(p2)
