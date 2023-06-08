
##############################################################
## Heade File for Topological Quantities from Quantum Circuits
##############################################################

include("../../QuantumSimulator_DEV/QSim.jl")

using DelimitedFiles
using LinearAlgebra
using Plots


##############
# Hamiltonians
##############

# Pauli matrices
sigma_x = [0. 1.; 1. 0.]
sigma_y = [0. -1.0im; 1.0im 0.]
sigma_z = [1. 0.; 0. -1.]


## 1D Hamiltonians ##

""" Momentum space bulk Hamiltonian of SSH (Rice-Mele) model. """
function SSH_model(kx, params)
   δ = params[1]
   m = params[2]
   tk = (1. + δ)*exp(1.0im*kx)
   return [-m 1.0+tk; 1.0+conj(tk) m]
end

""" Momentum space bulk Hamiltonian of QWZ model. """
function QWZ_model(kx, params)
   ky = params[1]
   u = params[2]
   return sin(kx)*sigma_x + sin(ky)*sigma_y + (u + cos(kx) + cos(ky))*sigma_z
end


""" Momentum space bulk Hamiltonian of 2D p-wave superconductor. """
function H_p_wave(kx, params)
    ky = params[1]
    Δ = params[2]
    t = params[3]
    μ = params[4]
    E_k = t*(cos(kx) + cos(ky)) + μ
    return Δ*(sin(ky)*sigma_x + sin(kx*sigma_y)) - E_k*sigma_z
end


## 2D Hamiltonians ##

""" Momentum space bulk Hamiltonian of QWZ model. """
function QWZ_model(kx, ky, params)
   u = params[1]
   return sin(kx)*sigma_x + sin(ky)*sigma_y + (u + cos(kx) + cos(ky))*sigma_z
end


""" Momentum space bulk Hamiltonian of the Haldane model. """
function Haldane_model(kx, ky, params)

    t1 = params[1]
    t2 = params[2]
    m = params[3]
    ϕ = params[4]
    a = params[5]

    # hardcoded, symmetrical basis
    a1_1 = sqrt(3)*a/2
    a1_2 = a/2
    a2_1 = sqrt(3)*a/2
    a2_2 = -a/2

    #a1_1 = params[5]
    #a1_2 = params[6]
    #a2_1 = params[7]
    #a2_2 = params[8]

    dx = t1*(cos(kx*a1_1+ky*a1_2) + cos(kx*a2_1+ky*a2_2) + 1)
    dy = t1*(sin(kx*a1_1+ky*a1_2) + sin(kx*a2_1+ky*a2_2))
    dz = m + 2*t2*sin(ϕ)*(sin(kx*a1_1+ky*a1_2) - sin(kx*a2_1+ky*a2_2) - sin(kx*(a1_1-a2_1)+ky*(a1_2-a2_2)))
    return dx*sigma_x + dy*sigma_y + dz*sigma_z
end


""" Momentum space bulk Hamiltonian of 2D p-wave superconductor. """
function H_p_wave(kx, ky, params)
    Δ = params[1]
    t = params[2]
    μ = params[3]
    E_k = t*(cos(kx) + cos(ky)) + μ
    return Δ*(sin(ky)*sigma_x + sin(kx*sigma_y)) - E_k*sigma_z
end


####################
# Exact Calculations
####################

# Functions to compute the Berry phase and Chern number of a given 1D or
# 2D Hamiltonian (numerically) exactly.


""" Find occupied states of a given Hamiltonian H for Berry phase calculation. """
function occupied_states(H)
   vals, wfs = eigen(H)
   occ_wfs = []
   for i in 1:length(vals)
      if Real(vals[i]) < 0
         push!(occ_wfs, wfs[:,i])
      end
   end
   return occ_wfs
end


""" Calculate the Berry phase by performing a closed loop [0, 2π] in
momentum space. Need to specify discretisation parameter N_ks. """
function Berry_phase(H, N_ks, params)

   # get momentum space [0, 2π]
   ks = LinRange(0, 1, N_ks)[1:N_ks-1]*2π

   # initialise
   wf0 = occupied_states(H(ks[1], params))
   wfold = copy(wf0)

   # loop through momentum space and update Berry phase
   local M = 1.
   for i in 2:N_ks-1
      wf = occupied_states(H(ks[i], params))
      M *= dot(wfold, wf)
      wfold = copy(wf)
   end
   M *= dot(wfold, wf0)
   Berry_phase = angle(M)

   return Berry_phase
end


""" Auxiliary function for Chern number calculation. Evaluates
the value of a single link in a plaquette, given a matrix of eigenvectors
for the whole BZ. """
function link(ulist, a1, b1, a2, b2, Nk)
    vec1 = ulist[mod(a1, 1:Nk), mod(b1, 1:Nk), :]
    vec2 = ulist[mod(a2, 1:Nk), mod(b2, 1:Nk), :]
    return dot(vec1, vec2)
end


""" Function to calculate the Chern number of a given momentum
space Hamiltonian by decomposing the Brillouin zone into Nk
plaquettes. """
function Chern_number(H, params, Nk, band_index=1)

    k_step = 2π/Nk
    nuc = size(H(0, 0, params))[1]

    # table of eigenvectors
    ulist = Complex.(zeros(Nk, Nk, nuc))
    for a in 1:Nk
        kx = -π + a*k_step
        for b in 1:Nk
            ky = -π + b*k_step
            ulist[a, b, :] = eigvecs(H(kx, ky, params))[:, band_index]
        end
    end

    # save Berry flux for each plaquette
    berryflux = Complex.(zeros(Nk, Nk))
    for a in 1:Nk
        for b in 1:Nk
            #println("Doing plaquette ($a, $b)")
            link1 = link(ulist, a, b, a+1, b, Nk)
            link2 = link(ulist, a+1, b, a+1, b+1, Nk)
            link3 = link(ulist, a+1, b+1, a, b+1, Nk)
            link4 = link(ulist, a, b+1, a, b, Nk)
            tmp = link1*link2*link3*link4
            bflux = log(tmp/abs(tmp))
            berryflux[a, b] = bflux
        end
    end

    # sum over all plaquettes to get Chern number
    return berryflux./(-1.0im), -real(sum(berryflux) / 2π*1.0im)
end


######################################
# Variational ground state preparation
######################################

# Unfinished functions to set up the ground state preparation (for
# a 2x2 Bloch Hamiltonian) variationally.


""" Function to calculate the expectation value of σx, σy or σz for a
1-qubit quantum circuit. """
function get_expectation_value(qc, op, N_meas)

    # rotate to right basis
    if op == "X"
        RyGate!(qc, [1], -π/2)

    elseif op == "Y"
        RxGate!(qc, [1], π/2)
    end
    sample_measurement(qc, [1], N_meas, plot=false, verbose=false)

    # evaluate expectation value
    meas_freqs = collect(values(qc.MeasurementResult))
    exp_value = (+1)*meas_freqs[1] + (-1)*meas_freqs[2]

    return exp_value
end


""" Function to get the estimated ground state energy of a Bloch Hamiltonian
of the form H = a*σx + b*σy + c*σz, with three rotation angles to be optimised. """
function energy_expectation_value_HBloch(params, a, b, c, N_meas; N_qubits=1, backend="MPS_ITensor", maxbond=2)

    # parameters to be optimised
    θ1 = params[1]
    θ2 = params[2]
    θ3 = params[3]

    # save expectation values and loop through circuit
    exp_vals = []
    for op in ["X", "Y", "Z"]

        # prepare quantum circuit
        local qc_var = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
        RxGate!(qc_var, [1], θ1)
        RzGate!(qc_var, [1], θ2)
        RxGate!(qc_var, [1], θ3)

        # measure each expectation value
        push!(exp_vals, get_expectation_value(qc, op, N_meas))

    end

    # construct energy estimate
    E_est = a*exp_vals[1] + b*exp_vals[2] + c*exp_vals[2]
    return E_est
end


""" Function to variationally minimise the ground state of the momentum
space QWZ Hamiltonian for each point in the Brillouin zone. """
function variational_ground_state_prep_QWZ(Nk, u; filled=true)

    # get Brillouin zone
    δk = 2π/Nk
    kxs = [-π + i*δk for i in 1:Nk]
    kys = [-π + i*δk for i in 1:Nk]

    optimised_params = zeros(Nk, Nk, 3)

    # loop through BZ
    for i in 1:Nk
        for j in 1:Nk

            # params for H at given point in BZ
            a = sin(kxs[i])
            b = sin(kys[j])
            c = u + cos(kxs[i]) + cos(kys[j])

            # do optimisation
            initial = [0., 0., 0.]
            result = optimize(params -> energy_expectation_value_HBloch(params, a, b, c, N_meas), initial, iterations=50000, allow_f_increases=true)


        end
    end

end


###########################################
# Quantum Circuit measuring the Berry Phase
###########################################

# Implementation of the adiabatic loop evolution to calculate the Berry
# phase (with cancellation of the dynamical phase). Quantum Circuits include
# both a Hadamard test circuit as well as a quantum phase estimation circuit.
# Still in progress to find correct results...

""" Function to construct the unitary similarity transforms for a given 1D
Bloch Hamiltonian at a given momentum in the BZ. """
function get_unitary_sim_transform(H, params)

    kx = 0
    eig_vecs = eigvecs(H(kx, params))
    filled_band = eig_vecs[:, 1]
    empty_band = eig_vecs[:, 2]

    U = Complex.(zeros(2, 2))
    U[:, 1] = filled_band[:]
    U[:, 2] = empty_band[ :]

    return U
end


""" Function to perform an adiabatic loop to calculate
the Berry phase. Starting point of the loop is given. Specifically
for Hadamard test. """
function adiabatic_loop_1D!(qc, H, params, Nk, dt)

    # evolution step (discretisation of BZ)
    δkx = 2π/(2*Nk)

    # positive time evolution
    for i in 1:Nk
        #kx_ev = -π + i*δkx
        kx_ev = 0 + i*δkx
        dU1 = exp(-1.0im*H(kx_ev, params)*dt)
        C_UGate!(qc, dU1, [1], [2])
    end

    # negative time evolution
    for i in 1:Nk
        #kx_ev = -π + Nk*δkx + i*δkx
        kx_ev = π + i*δkx
        dU2 = exp(1.0im*H(kx_ev, params)*dt)
        C_UGate!(qc, dU2, [1], [2])
    end
end


""" Function to perform an adiabatic loop to calculate
the Berry phase. Starting point of the loop is given.  More general
for QPE type algorithm. """
function adiabatic_loop_1D!(qc, H, params, control, pos, Nk, dt)

    # evolution step (discretisation of BZ)
    δkx = 2π/(2*Nk)

    # positive time evolution
    for i in 1:Nk
        #kx_ev = -π + i*δkx
        kx_ev = 0 + i*δkx
        dU1 = exp(-1.0im*H(kx_ev, params)*dt)
        C_UGate!(qc, dU1, control, pos)
    end

    # negative time evolution
    for i in 1:Nk
        #kx_ev = -π + Nk*δkx + i*δkx
        kx_ev = π + i*δkx
        dU2 = exp(1.0im*H(kx_ev, params)*dt)
        C_UGate!(qc, dU2, control, pos)
    end
end



""" Function to set up a Hadamard test circuit for the Berry phase of
a 1D Hamiltonian. Can select whether a "cos" or a "sin" circuit should
be calculated. """
function Berry_phase_Hadamard_test(H, params, U_init, N_meas, Nk, dt; cos_qc=true, N_qubits=2, backend="MPS_ITensor", maxbond=10)

    # quantum circuit for cos
    if cos_qc
        qc_cos = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
        Hadamard!(qc_cos, [1])
        UGate!(qc_cos, U_init, [2])
        adiabatic_loop_1D!(qc_cos, H, params, Nk, dt)
        Hadamard!(qc_cos, [1])
        sample_measurement(qc_cos, [1], N_meas, plot=false, verbose=false, eps=0.000001)
        meas_freqs_cos = collect(values(qc_cos.MeasurementResult))
        #return meas_freqs_cos[1] # returning correct result?
        println("meas_freqs_cos ", meas_freqs_cos)
        return meas_freqs_cos # returning correct result?

    # quantum circuit for sine
    else
        qc_sin = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
        Hadamard!(qc_sin, [1])
        SGate!(qc_sin, [1]) # extra S-gate
        UGate!(qc_sin, U_init, [2])
        adiabatic_loop_1D!(qc_sin, H, params, Nk, dt)
        Hadamard!(qc_sin, [1])
        sample_measurement(qc_sin, [1], N_meas, plot=false, verbose=false, eps=0.000001)
        meas_freqs_sin = collect(values(qc_sin.MeasurementResult))
        #return meas_freqs_sin[1] # returning correct result?
        println("meas_freqs_sin ", meas_freqs_sin)
        return meas_freqs_sin
    end
end


""" Function to measure the Berry phase in 1D via an adiabatic loop. Constructs
two Hadamard test circuits from which the Berry is calculated. """
function measure_Berry_phase_1D(H, params, N_meas, Nk, dt; N_qubits=2, backend="MPS_ITensor", maxbond=10)

    # unitary transformation to initial state
    U_init = get_unitary_sim_transform(H, params)

    # get measurement results
    P_cos = Berry_phase_Hadamard_test(H, params, U_init, N_meas, Nk, dt, cos_qc=true, N_qubits=N_qubits, backend=backend, maxbond=maxbond)
    P_sin = Berry_phase_Hadamard_test(H, params, U_init, N_meas, Nk, dt, cos_qc=false, N_qubits=N_qubits, backend=backend, maxbond=maxbond)

    # is the order of the probabilities here right?
    #P_0_cos = P_cos[1]
    #P_1_cos = P_cos[2]
    #P_0_sin = P_sin[1]
    #P_1_sin = P_sin[2]

    P_0_cos = P_cos[2]
    P_1_cos = P_cos[1]
    P_0_sin = P_sin[1]
    P_1_sin = P_sin[2]

    # get sin and cos of theta
    cos_θ = 2*P_0_cos - 1
    sin_θ = 1 - 2*P_0_sin

    # recover Berry phase based on interval where inverse function is defined
    if cos_θ <= 0 && sin_θ <= 0 # -π < θ < -π/2
        println("region 1")
        return -π + abs(asin(1 - 2*P_0_sin))
    elseif cos_θ >= 0 && sin_θ <= 0 # -π/2 < θ < 0
        println("region 2")
        return asin(1 - 2*P_0_sin)
    elseif cos_θ >= 0 && sin_θ >= 0 # 0 < θ < π/2
        println("region 3")
        return acos(2*P_0_cos - 1)
    elseif cos_θ <= 0 && sin_θ >= 0 # π/2 < θ < π
        println("region 4")
        return acos(2*P_0_cos - 1)
    end
end


""" Function to calculate the Berry phase via a quantum phase estimation
algorithm. """
function Berry_phase_QPE(H, params, Nk, dt, N_qubits, n_prec, N_meas; backend="MPS_ITensor", maxbond=20)

    # get quantum circuit
    qc = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)

    # initialise phase estimation register
    Hadamard!(qc, [i for i in 1:N_qubits-1])

    # prepare initial state (kx = 0)
    U_init = get_unitary_sim_transform(H, params)
    UGate!(qc, U_init, [N_qubits])

    # add extra phase to |1⟩ to map Berry phase into correct range
    local m = 1
    for i in (N_qubits-1):-1:(1)
        for loop in 1:m
            SGate!(qc, [i])
            SGate!(qc, [i])
        end
        m = 2*m
    end

    # do controlled rotations
    local n = 1
    for i in (N_qubits-1):-1:(1)
        println("Doing i = $i and n = $n")
        for loop in 1:n
            adiabatic_loop_1D!(qc, H, params, [i], [N_qubits], Nk, dt)
        end
        n = 2*n
    end

    # do inverse QFT (switch off representation)
    invQFT!(qc, 1, N_qubits-1)

    # measurement
    sample_measurement(qc, [i for i in 1:n_prec], N_meas, plot=false, verbose=true, eps=0.000001)
    phases, probs_max = QPE_get_phase(qc, 10)

    #draw(qc)

    return phases, probs_max
    #return phases[1]
end


""" TEST CODE """
function QPE_test(H, params, Nk, dt, N_qubits, n_prec, N_meas; backend="MPS_ITensor", maxbond=20)

    # get quantum circuit
    qc = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)

    # initialise phase estimation register
    Hadamard!(qc, [i for i in 1:N_qubits-1])

    # prepare initial state (kx = 0)
    U_init = get_unitary_sim_transform(H, params)
    UGate!(qc, U_init, [N_qubits])

    # do controlled rotations
    local n = 1
    for i in (N_qubits-1):-1:(1)
        println("Doing i = $i and n = $n")
        for loop in 1:n
            #adiabatic_loop_1D!(qc, H, params, [i], [N_qubits], Nk, dt)
            e_iH = exp(1.0im*H(0, params))
            C_UGate!(qc, e_iH, [i], [N_qubits])
        end
        n = 2*n
    end

    # do inverse QFT (switch off representation)
    invQFT!(qc, 1, N_qubits-1)

    # measurement
    sample_measurement(qc, [i for i in 1:n_prec], N_meas, plot=true, verbose=true, eps=0.000001)
    phases, probs_max = QPE_get_phase(qc, 3)

    #draw(qc)

    return phases, probs_max
end


##########################################
# Quantum Circuit measuring distinct links
##########################################

# Quantum circuits which measure the Chern number by decomposing the
# Brillouin zone into plaquettes and calculating the expectation value
# of each link separately. No adiabatic loop involved here.


""" Function to construct a list of unitary similarity transforms for a given
Bloch Hamiltonian in a given 2D Brillouin zone. """
function get_unitary_sim_transform(H, params, Nk; filled=true)

    # get Brillouin zone
    δk = 2π/Nk
    kxs = [-π + i*δk for i in 1:Nk]
    kys = [-π + i*δk for i in 1:Nk]

    # table of eigenvectors
    filled_band_list = Complex.(zeros(Nk, Nk, 2))
    empty_band_list = Complex.(zeros(Nk, Nk, 2))
    for a in 1:Nk
        for b in 1:Nk
            eig_vecs = eigvecs(H(kxs[a], kys[b], params))
            filled_band_list[a, b, :] = eig_vecs[:, 1]
            empty_band_list[a, b, :] = eig_vecs[:, 2]
        end
    end

    # construct unitary similarity transformations over BZ
    U_list = Complex.(zeros(Nk, Nk, 2, 2))
    for a in 1:Nk
        for b in 1:Nk
            U = Complex.(zeros(2, 2))
            U[:, 1] = filled_band_list[a, b, :]
            U[:, 2] = empty_band_list[a, b, :]
            U_list[a, b, :, :] = U
        end
    end

    return U_list
end


""" Function to peform Chern number computation on a 2 qubit system: sets
up a Hadamard test circuit to measure either the real or imaginary value
of a single link in a plaquette. """
function Hadamard_test_Chern_number!(qc, U, U_ev, params; filled=true)

    # set up, transform to initial state in BZ via U
    Hadamard!(qc, [1])
    UGate!(qc, U, [2])

    # controlled evolution to kx_ev, ky_ev state
    C_UGate!(qc, U', [1], [2])
    C_UGate!(qc, U_ev, [1], [2])

end


""" Function to evaluate the (real or imaginary part of the) expectation
value of the operator U in the Chern number circuit. """
function U_expectation_value!(qc, N_meas; Re=true)

    # do corresponding measurement
    if Re # measure in x basis
        RyGate!(qc, [1], -π/2)
        #Hadamard!(qc, [1])
    else # Im, measure in y basis
        RxGate!(qc, [1], π/2)
        #S_dagGate!(qc, [1])
        #Hadamard!(qc, [1])
    end
    sample_measurement(qc, [1], N_meas, plot=false, verbose=false)

    # evaluate expectation value
    meas_freqs = collect(values(qc.MeasurementResult))
    exp_value = (+1)*meas_freqs[1] + (-1)*meas_freqs[2]

    return exp_value
end


""" Function to measure the Chern number of a given Hamiltonian. Decomposes the
Brillouin zone into plaquettes, calculates each link via a 2-qubit quantum circuit
and obtains the Berry flux for each plaquette by multuplying the links. """
function measure_Chern_number(H, params, Nk, N_meas; N_qubits=2, backend="MPS_ITensor", maxbond=10)

    # get Brillouin zone
    #δ = 2π/Nk
    #kxs = [-π + i*δk for i in 1:Nk]
    #kys = [-π + i*δk for i in 1:Nk]

    # save Berry fluxes
    Berry_flux = Complex.(zeros(Nk, Nk))

    # loop through BZ to find unitary similarity transforms
    U_list = get_unitary_sim_transform(H, params, Nk)

    # loop through BZ
    for i in 1:Nk # kx
        for j in 1:Nk # ky
            println("Doing plaquette ($(i), $(j))")

            ## first link ##
            U = U_list[mod(i, 1:Nk), mod(j, 1:Nk), :, :]
            U_ev = U_list[mod(i+1, 1:Nk), mod(j, 1:Nk), :, :]

            qc_re = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
            Hadamard_test_Chern_number!(qc_re, U, U_ev, params)
            U_re = U_expectation_value!(qc_re, N_meas, Re=true)

            qc_im = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
            Hadamard_test_Chern_number!(qc_im, U, U_ev, params)
            U_im = U_expectation_value!(qc_im, N_meas, Re=false)

            # reconstruct expectation value of link
            U_k1 = (U_re + 1.0im*U_im)/abs(U_re + 1.0im*U_im)


            ## second link ##
            U = U_list[mod(i+1, 1:Nk), mod(j, 1:Nk), :, :]
            U_ev = U_list[mod(i+1, 1:Nk), mod(j+1, 1:Nk), :, :]

            qc_re = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
            Hadamard_test_Chern_number!(qc_re, U, U_ev, params)
            U_re = U_expectation_value!(qc_re, N_meas, Re=true)

            qc_im = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
            Hadamard_test_Chern_number!(qc_im, U, U_ev, params)
            U_im = U_expectation_value!(qc_im, N_meas, Re=false)

            # reconstruct expectation value of link
            U_k2 = (U_re + 1.0im*U_im)/abs(U_re + 1.0im*U_im)


            ## third link ##
            U = U_list[mod(i+1, 1:Nk), mod(j+1, 1:Nk), :, :]
            U_ev = U_list[mod(i, 1:Nk), mod(j+1, 1:Nk), :, :]

            qc_re = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
            Hadamard_test_Chern_number!(qc_re, U, U_ev, params)
            U_re = U_expectation_value!(qc_re, N_meas, Re=true)

            qc_im = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
            Hadamard_test_Chern_number!(qc_im, U, U_ev, params)
            U_im = U_expectation_value!(qc_im, N_meas, Re=false)

            # reconstruct expectation value of link
            U_k3 = (U_re + 1.0im*U_im)/abs(U_re + 1.0im*U_im)


            ## fourth link ##
            U = U_list[mod(i, 1:Nk), mod(j+1, 1:Nk), :, :]
            U_ev = U_list[mod(i, 1:Nk), mod(j, 1:Nk), :, :]

            qc_re = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
            Hadamard_test_Chern_number!(qc_re, U, U_ev, params)
            U_re = U_expectation_value!(qc_re, N_meas, Re=true)

            qc_im = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
            Hadamard_test_Chern_number!(qc_im, U, U_ev, params)
            U_im = U_expectation_value!(qc_im, N_meas, Re=false)

            # reconstruct expectation value of link
            U_k4 = (U_re + 1.0im*U_im)/abs(U_re + 1.0im*U_im)

            # evaluate Berry flux -> should be imaginary!
            Berry_flux[i, j] = log(U_k1*U_k2*U_k3*U_k4)

        end
    end

    # evaluate Chern number
    return real(sum(Berry_flux)/(2π*1.0im))
end



#################################################################
# Quantum Circuit performing adiabatic evolution around plaquette
#################################################################

# Quntum circuits which calculate the Chern number by decomposing the
# Brillouin zone into plaquettes and evaluating the Berry phase around
# each plaquette via an adiabatic (double) loop. More in the spirit of
# what can actually be done on quantum computers.


""" Function to perform an adiabatic double loop around a plaquette to calculate
the Berry flux. Starting point of the plaquette is given. On the first loop,
the time evolution is positive and on the second loop it is negative,
such that the dynamical phase is exactly cancelled and only the Berry
phase remains (read out via Hadamard test or QPE)."""
function adiabatic_plaquette_loop!(qc, H, kx, ky, params, N_link, δk, dt)

    # evolution step along links of plaquette
    δk_link = δk/N_link


    ## first loop, positive time evolution ##

    # kx, ky -> kx + δk, ky
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + i*δk_link
        ky_ev = ky
        dU1 = exp(-1.0im*H(kx_ev, ky_ev, params)*dt) # positive time evolution
        C_UGate!(qc, dU1, [1], [2])
    end

    # kx + δk, ky -> kx + δk, ky + δk
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + N_link*δk_link
        ky_ev = ky + i*δk_link
        dU2 = exp(-1.0im*H(kx_ev, ky_ev, params)*dt) # positive time evolution
        C_UGate!(qc, dU2, [1], [2])
    end

    # kx + δk, ky + δk -> kx, ky + δk
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + N_link*δk_link - i*δk_link
        ky_ev = ky + N_link*δk_link
        dU3 = exp(-1.0im*H(kx_ev, ky_ev, params)*dt) # positive time evolution
        C_UGate!(qc, dU3, [1], [2])
    end

    # kx, ky + δk -> kx, ky
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx
        ky_ev = ky + N_link*δk_link - i*δk_link
        dU4 = exp(-1.0im*H(kx_ev, ky_ev, params)*dt) # positive time evolution
        C_UGate!(qc, dU4, [1], [2])
    end


    ## second loop, negative time evolution ##

    # kx, ky -> kx + δk, ky
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + i*δk_link
        ky_ev = ky
        dU1 = exp(1.0im*H(kx_ev, ky_ev, params)*dt) # negative time evolution
        C_UGate!(qc, dU1, [1], [2])
    end

    # kx + δk, ky -> kx + δk, ky + δk
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + N_link*δk_link
        ky_ev = ky + i*δk_link
        dU2 = exp(1.0im*H(kx_ev, ky_ev, params)*dt) # negative time evolution
        C_UGate!(qc, dU2, [1], [2])
    end

    # kx + δk, ky + δk -> kx, ky + δk
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + N_link*δk_link - i*δk_link
        ky_ev = ky + N_link*δk_link
        dU3 = exp(1.0im*H(kx_ev, ky_ev, params)*dt) # negative time evolution
        C_UGate!(qc, dU3, [1], [2])
    end

    # kx, ky + δk -> kx, ky
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx
        ky_ev = ky + N_link*δk_link - i*δk_link
        dU4 = exp(1.0im*H(kx_ev, ky_ev, params)*dt) # negative time evolution
        C_UGate!(qc, dU4, [1], [2])
    end
end



""" Function to perform an adiabatic double loop around a plaquette to calculate
the Berry flux. Starting point of the plaquette is given. On the first loop,
the time evolution is positive and on the second loop it is negative,
such that the dynamical phase is exactly cancelled and only the Berry
phase remains (read out via Hadamard test or QPE)."""
function adiabatic_plaquette_loop!(qc, H, kx, ky, params, control, pos, N_link, δk, dt)

    # evolution step along links of plaquette
    δk_link = δk/N_link


    ## first loop, positive time evolution ##

    # kx, ky -> kx + δk, ky
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + i*δk_link
        ky_ev = ky
        dU1 = exp(-1.0im*H(kx_ev, ky_ev, params)*dt) # positive time evolution
        C_UGate!(qc, dU1, control, pos)
    end

    # kx + δk, ky -> kx + δk, ky + δk
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + N_link*δk_link
        ky_ev = ky + i*δk_link
        dU2 = exp(-1.0im*H(kx_ev, ky_ev, params)*dt) # positive time evolution
        C_UGate!(qc, dU2, control, pos)
    end

    # kx + δk, ky + δk -> kx, ky + δk
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + N_link*δk_link - i*δk_link
        ky_ev = ky + N_link*δk_link
        dU3 = exp(-1.0im*H(kx_ev, ky_ev, params)*dt) # positive time evolution
        C_UGate!(qc, dU3, control, pos)
    end

    # kx, ky + δk -> kx, ky
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx
        ky_ev = ky + N_link*δk_link - i*δk_link
        dU4 = exp(-1.0im*H(kx_ev, ky_ev, params)*dt) # positive time evolution
        C_UGate!(qc, dU4, control, pos)
    end


    ## second loop, negative time evolution ##

    # kx, ky -> kx + δk, ky
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + i*δk_link
        ky_ev = ky
        dU1 = exp(1.0im*H(kx_ev, ky_ev, params)*dt) # negative time evolution
        C_UGate!(qc, dU1, control, pos)
    end

    # kx + δk, ky -> kx + δk, ky + δk
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + N_link*δk_link
        ky_ev = ky + i*δk_link
        dU2 = exp(1.0im*H(kx_ev, ky_ev, params)*dt) # negative time evolution
        C_UGate!(qc, dU2, control, pos)
    end

    # kx + δk, ky + δk -> kx, ky + δk
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx + N_link*δk_link - i*δk_link
        ky_ev = ky + N_link*δk_link
        dU3 = exp(1.0im*H(kx_ev, ky_ev, params)*dt) # negative time evolution
        C_UGate!(qc, dU3, control, pos)
    end

    # kx, ky + δk -> kx, ky
    for i in 1:N_link # count every point around plaquette only once
        kx_ev = kx
        ky_ev = ky + N_link*δk_link - i*δk_link
        dU4 = exp(1.0im*H(kx_ev, ky_ev, params)*dt) # negative time evolution
        C_UGate!(qc, dU4, control, pos)
    end
end


""" Function to perform a Hadamard test to evaluate the Berry flux for a
given plaquette via an adiabatic double loop around the plaquette. """
function Hadamard_test_adiabatic_loop(H, params, U, kx, ky, N_link, δk, dt, N_meas; N_qubits=2, backend="MPS_ITensor", maxbond=10)

    # initialise quantum circuit, set to eigenstate in plaquette
    qc_plaq = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)
    Hadamard!(qc_plaq, [1])
    SGate!(qc_plaq, [1]) # extra S-gate to get sine in the end
    #U = U_list[mod(i, 1:Nk), mod(j, 1:Nk), :, :]
    UGate!(qc_plaq, U, [2])

    # peform adiabatic loop
    adiabatic_plaquette_loop!(qc_plaq, H, kx, ky, params, N_link, δk, dt)

    # wrap up and measure
    Hadamard!(qc_plaq, [1])
    sample_measurement(qc_plaq, [1], N_meas, plot=false, verbose=true, eps=0.000001)

    # evaluate Berry flux
    meas_freqs_sin = collect(values(qc_plaq.MeasurementResult))

    P_0_sin = meas_freqs_sin[1]
    #P_0_sin = meas_freqs_sin[2]
    println("sin(2θ) = ", 1 - 2*P_0_sin)
    return 0.5*asin(1 - 2*P_0_sin)
end


""" Function implementing the Berry phase calculation for a given plaquette
via an adiabatic loop with a QPE.  """
function adiabatic_loop_QPE(H, params, kx, ky, dt, N_link, δk, N_qubits, n_prec, N_meas; backend="MPS_ITensor", maxbond=20)

    # get quantum circuit
    qc = initialise_qcircuit(N_qubits, backend, maxdim=maxbond)

    # initialise phase estimation register
    Hadamard!(qc, [i for i in 1:N_qubits-1])

    # add extra phase to |1⟩ to map Berry phase into correct range
    local m = 1
    for i in (N_qubits-1):-1:(1)
        for loop in 1:m
            SGate!(qc, [i])
            SGate!(qc, [i])
        end
        m = 2*m
    end

    # do controlled rotations
    local n = 1
    for i in (N_qubits-1):-1:(1)
        println("Doing i = $i and n = $n")
        for loop in 1:n
            adiabatic_plaquette_loop!(qc, H, kx, ky, params, [i], [N_qubits], N_link, δk, dt)
        end
        n = 2*n
    end

    # do inverse QFT (switch off representation)
    invQFT!(qc, 1, N_qubits-1)

    # measurement
    sample_measurement(qc, [i for i in 1:n_prec], N_meas, plot=false, verbose=false, eps=0.000001)
    phases, probs_max = QPE_get_phase(qc, 3)

    # return Berry phase correctly
    #if phases[1] < π
    #     return phases[1]
    #else
    #     return phases[1] - 2π
    #end

    return phases[1]-π
end


""" Function to measure the Chern number of a given Hamiltonian by decomposing
the Brillouin zone into plaquettes and performing an Hadamard test via a
double adiabatic loop for each plaquette to find the corresponding Berry flux. """
function measure_Chern_number_adiabatic_loop(H, params, N_link, Nk, dt, N_meas)

    # Brillouin zone spacing
    δk = 2π/Nk

    # save Berry fluxes
    Berry_flux = Complex.(zeros(Nk, Nk))

    # loop through BZ to find unitary similarity transforms
    U_list = get_unitary_sim_transform(H, params, Nk)

    # loop through BZ
    for i in 1:Nk # kx
        for j in 1:Nk # ky
            println("Doing plaquette ($(i), $(j))")
            U = U_list[mod(i, 1:Nk), mod(j, 1:Nk), :, :]
            kx = -π + i*δk
            ky = -π + j*δk
            Berry_flux[i, j] = Hadamard_test_adiabatic_loop(H, params, U, kx, ky, N_link, δk, dt, N_meas)
        end
    end
    println("Berry flux: ", Berry_flux)

    # evaluate Chern number
    #return real(sum(Berry_flux)/(2π*1.0im))
    return Berry_flux, real(sum(Berry_flux)/(2π))
end


""" Function to measure the Chern number of a given Hamiltonian by decomposing
the Brillouin zone into plaquettes and performing a QPE  via a double adiabatic
loop for each plaquette to find the corresponding Berry flux. """
function measure_Chern_number_adiabatic_loop_QPE(H, params, N_link, Nk, dt, N_qubits, n_prec, N_meas; backend="MPS_ITensor", maxdim=20)

    # Brillouin zone spacing
    δk = 2π/Nk

    # save Berry fluxes
    Berry_flux = Complex.(zeros(Nk, Nk))

    # loop through BZ to find unitary similarity transforms
    U_list = get_unitary_sim_transform(H, params, Nk)

    # loop through BZ
    for i in 1:Nk # kx
        for j in 1:Nk # ky
            println("Doing plaquette ($(i), $(j))")
            U = U_list[mod(i, 1:Nk), mod(j, 1:Nk), :, :]
            kx = -π + i*δk
            ky = -π + j*δk
            Berry_flux[i, j] = adiabatic_loop_QPE(H, params, kx, ky, dt, N_link, δk, N_qubits, n_prec, N_meas; backend=backend, maxbond=maxdim)
        end
    end

    # evaluate Chern number
    return real(sum(Berry_flux)/(2π*1.0im))
end


###########
# old stuff
###########


function measure_Chern_number_plaquette_comparison(H, params, N_link, Nk, dt, N_meas)

    # Brillouin zone spacing
    δk = 2π/Nk

    # save Berry fluxes
    #Berry_flux_ex = Complex.(zeros(Nk, Nk))
    #Berry_flux = Complex.(zeros(Nk, Nk))
    Berry_flux_ex = zeros(Nk, Nk)
    Berry_flux = zeros(Nk, Nk)

    # loop through BZ to find unitary similarity transforms
    U_list = get_unitary_sim_transform(H, params, Nk)

    # table of eigenvectors
    nuc = size(H(0, 0, params))[1]
    band_index=1
    ulist = Complex.(zeros(Nk, Nk, nuc))
    for a in 1:Nk
        kx = -π + a*δk
        for b in 1:Nk
            ky = -π + b*δk
            ulist[a, b, :] = eigvecs(H(kx, ky, params))[:, band_index]
        end
    end

    # loop through BZ
    for i in 1:Nk # kx
        for j in 1:Nk # ky
            println("Doing plaquette ($(i), $(j))")

            # adiabatic loop
            U = U_list[mod(i, 1:Nk), mod(j, 1:Nk), :, :]
            kx = -π + i*δk
            ky = -π + j*δk
            Berry_flux[i, j] = real(Hadamard_test_adiabatic_loop(H, params, U, kx, ky, N_link, δk, dt, N_meas))

            # exact
            link1 = link(ulist, i, j, i+1, j, Nk)
            link2 = link(ulist, i+1, j, i+1, j+1, Nk)
            link3 = link(ulist, i+1, j+1, i, j+1, Nk)
            link4 = link(ulist, i, j+1, i, j, Nk)
            tmp = link1*link2*link3*link4
            Berry_flux_ex[i, j] = real(log(tmp/abs(tmp))/(-1.0im))

            println("Flux exact: $(Berry_flux_ex[i, j])")
            println("Flux adiabatic loop: $(Berry_flux[i, j])")
            println("equal? ", Berry_flux_ex[i, j]==Berry_flux[i, j])
        end
    end

    println("Berry flux: ", Berry_flux)

    # evaluate Chern number
    #return real(sum(Berry_flux)/(2π*1.0im))
    return Berry_flux, Berry_flux_ex, real(sum(Berry_flux)/(2π)), real(sum(Berry_flux_ex)/2π)
end
