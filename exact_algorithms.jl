
#######################################################
## Exact (analytical) evaluations of quantum algorithms
#######################################################


##############
# Hamiltonians
##############


""" Function to obtain an N-site tensor product of the operator O,
as E ⊗ ... ⊗ O ⊗ E ⊗ ... ⊗ E. """
function on_site_operator(O, N, i)

    s = Index(2, "QCircuit")
    E = array(op("E", s))

    # initialise
    if i == 1
        on_site = O
    else
        on_site = E
    end

    for j in 2:N
        if j == i
            on_site = kron(on_site, O)
        else
            on_site = kron(on_site, E)
        end
    end
    return on_site
end



""" Function to obtain the matrix representing the N-site operator
E ⊗ ... ⊗ O1 ⊗ O2 ⊗ E ⊗ ... ⊗ E., where O1 and O2 are in positions i and
i+1. """
function nearest_neighbor_coupling(O1, O2, N, i)

    if i < N
        op1 = on_site_operator(O1, N, i)
        op2 = on_site_operator(O2, N, i+1)

    else # i == N, PBC case
        op1 = on_site_operator(O1, N, N)
        op2 = on_site_operator(O2, N, 1)
    end

    return op1*op2
end


""" Function to obtain the exact matrix representation of an N-site
(anisotropic, transverse field) Heisenberg model. """
function Heisenberg_model(N, Jx=1, Jy=1, Jz=1, h=0, mag="FM", PBC=true)

    s = Index(2, "QCircuit")
    X = array(op("X", s))
    Y = array(op("Y", s))
    Z = array(op("Z", s))

    # initialise Hamiltonian
    Hxx = nearest_neighbor_coupling(X, X, N, 1)
    Hyy = nearest_neighbor_coupling(Y, Y, N, 1)
    Hzz = nearest_neighbor_coupling(Z, Z, N, 1)
    Hz = on_site_operator(Z, N, 1)

    # construct terms in Hamiltonian
    for i in 2:N
        if i < N
            Hxx += nearest_neighbor_coupling(X, X, N, i)
            Hyy += nearest_neighbor_coupling(Y, Y, N, i)
            Hzz += nearest_neighbor_coupling(Z, Z, N, i)
        elseif i == N && PBC
            Hxx += nearest_neighbor_coupling(X, X, N, i)
            Hyy += nearest_neighbor_coupling(Y, Y, N, i)
            Hzz += nearest_neighbor_coupling(Z, Z, N, i)
        end
        Hz += on_site_operator(Z, N, 2)
    end

    # construct and retun full Hamiltonian
    if mag == "FM"
        return -0.5*(Jx*Hxx + Jy*Hyy + Jz*Hzz) + h*Hz
    else # AFM
        return 0.5*(Jx*Hxx + Jy*Hyy + Jz*Hzz) + h*Hz
    end
end



###########################
# Quantum Fourier transform
###########################


""" Function to calculate the coefficients of a Fourier-transformed
state directly. """
function QFT_ex(t, initial_state)

    coeffs = Complex.(zeros(2^t))
    ω = exp(1.0im*2*π/(2.0^t))

    for j in 0:2^t-1
        b_j = Complex(0.0)
        for k in 0:(2^t-1)
            b_j += initial_state[k+1] * ω^(j*k)
        end
        coeffs[j+1] = b_j/√(2^t)
    end
    return coeffs
end


""" Function to calculate the coefficients of an inverse Fourier-transformed
state directly. """
function invQFT_ex(t, initial_state)

    coeffs = Complex.(zeros(2^t))
    ω = exp(-1.0im*2*π/(2.0^t))

    for j in 0:2^t-1
        b_j = Complex(0.0)
        for k in 0:(2^t-1)
            b_j += initial_state[k+1] * ω^(j*k)
        end
        coeffs[j+1] = b_j/√(2^t)
    end
    return coeffs
end



##########################
# Quantum phase estimation
##########################


""" Function to evaluate the coefficients of the quantum phase estimation
exactly for the two qubit case. """
function phase_estimation_register_exact(θ, t)

    coeffs = [1, exp(2π*1.0im * θ * 2^(t-1))]
    for i in t-2:-1:0
        coeffs = kron(coeffs, [1, exp(2π*1.0im * θ * 2^i)])
    end

    return invQFT_ex(t, 1/2^(t/2) * coeffs)
end


""" Function to combine the exactly simulated register of a QPE with
the eigenstate that is input in the algorithm. """
function phase_estimation_state_exact(register, eigenstate, t)
    return kron(register, eigenstate)
end
