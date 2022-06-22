
############################################
## Hilbert space for N-qubit quantum circuit
############################################

using ITensors


ITensors.space(::SiteType"QCircuit") = 2
#ITensors.space(::SiteType"QC") = 2

"""
function ITensors.space(
  ::SiteType"QCircuit";
  conserve_qns=false,
  conserve_parity=conserve_qns,
  conserve_number=false,
  qnname_parity="Parity",
  qnname_number="Number",
)
  if conserve_number && conserve_parity
    return [
      QN((qnname_number, 0), (qnname_parity, 0, 2)) => 1,
      QN((qnname_number, 1), (qnname_parity, 1, 2)) => 1,
    ]
  elseif conserve_number
    return [QN(qnname_number, 0) => 1, QN(qnname_number, 1) => 1]
  elseif conserve_parity
    return [QN(qnname_parity, 0, 2) => 1, QN(qnname_parity, 1, 2) => 1]
  end
  return 2
end
"""

# what the hell is wrong with that?
#ITensors.val(::ValName"0", ::SiteType"QCircuit") = 0
#ITensors.val(::ValName"1", ::SiteType"QCircuit") = 1

ITensors.state(::StateName"0", ::SiteType"QCircuit") = [1.0, 0.0]
ITensors.state(::StateName"1", ::SiteType"QCircuit") = [0.0, 1.0]


###################
# single-site gates
###################


function ITensors.op!(Op::ITensor, ::OpName"E",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>1] = Complex(1.)
    Op[s'=>2, s=>2] = Complex(1.)

end


function ITensors.op!(Op::ITensor, ::OpName"Proj00",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>1] = Complex(1.)

end


function ITensors.op!(Op::ITensor, ::OpName"Proj11",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>2, s=>2] = Complex(1.)

end


# |0><1|: takes spin down and turns into spin up
function ITensors.op!(Op::ITensor, ::OpName"S+",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>2] = Complex(1.)

end


# |1><0|: takes spin up and turns into spin down
function ITensors.op!(Op::ITensor, ::OpName"S-",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>2, s=>1] = Complex(1.)

end


function ITensors.op!(Op::ITensor, ::OpName"X",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>2] = Complex(1.)
    Op[s'=>2, s=>1] = Complex(1.)

end


function ITensors.op!(Op::ITensor, ::OpName"Y",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>2] = -1.0im
    Op[s'=>2, s=>1] = 1.0im

end


function ITensors.op!(Op::ITensor, ::OpName"Z",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>1] = Complex(1.)
    Op[s'=>2, s=>2] = Complex(-1.)

end


function ITensors.op!(Op::ITensor, ::OpName"H",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>1] = Complex(1.0 * 1/sqrt(2))
    Op[s'=>1, s=>2] = Complex(1.0 * 1/sqrt(2))
    Op[s'=>2, s=>1] = Complex(1.0 * 1/sqrt(2))
    Op[s'=>2, s=>2] = Complex(-1.0 * 1/sqrt(2))

end


function ITensors.op!(Op::ITensor, ::OpName"S",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>1] = Complex(1.0)
    Op[s'=>2, s=>2] = 1.0im

end


function ITensors.op!(Op::ITensor, ::OpName"T",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>1] = Complex(1.0)
    Op[s'=>2, s=>2] = exp(1.0im*pi/4)

end


function ITensors.op!(Op::ITensor, ::OpName"√X",
    ::SiteType"QCircuit", s::Index)

    Op[s'=>1, s=>1] = 1. + 1.0im
    Op[s'=>1, s=>2] = 1.0 - 1.0im
    Op[s'=>2, s=>1] = 1.0 - 1.0im
    Op[s'=>2, s=>2] = -1.0 + 1.0im

end


function ITensors.op!(Op::ITensor, ::OpName"Rx",
    ::SiteType"QCircuit", s::Index; θ::Number)

    Op[s'=>1, s=>1] = Complex(cos(θ/2))
    Op[s'=>1, s=>2] = -im*sin(θ/2)
    Op[s'=>2, s=>1] = -im*sin(θ/2)
    Op[s'=>2, s=>2] = Complex(cos(θ/2))

end


function ITensors.op!(Op::ITensor, ::OpName"Ry",
    ::SiteType"QCircuit", s::Index; θ::Number)

    Op[s'=>1, s=>1] = Complex(cos(θ/2))
    Op[s'=>1, s=>2] = -sin(θ/2)
    Op[s'=>2, s=>1] = sin(θ/2)
    Op[s'=>2, s=>2] = Complex(cos(θ/2))

end


function ITensors.op!(Op::ITensor, ::OpName"Rz",
    ::SiteType"QCircuit", s::Index; θ::Number)

    Op[s'=>1, s=>1] = Complex(1.0)
    Op[s'=>2, s=>2] = exp(im*θ)

end


# for quantum Fourier transform
function ITensors.op!(Op::ITensor, ::OpName"Rn",
    ::SiteType"QCircuit", s::Index; n::Number)

    Op[s'=>1, s=>1] = Complex(1.0)
    Op[s'=>2, s=>2] = exp(sign(n)*2π*1.0im/2^abs(n))

end


function ITensors.op!(Op::ITensor, ::OpName"P",
    ::SiteType"QCircuit", s::Index; θ::Number)

    Op[s'=>1, s=>1] = exp(1.0im * θ)
    Op[s'=>2, s=>2] = exp(1.0im * θ)

end


function ITensors.op!(Op::ITensor, ::OpName"U",
    ::SiteType"QCircuit", s::Index; α::Number, β::Number, γ::Number, δ::Number)

    Op[s'=>1, s=>1] = α
    Op[s'=>1, s=>2] = β
    Op[s'=>2, s=>1] = γ
    Op[s'=>2, s=>2] = δ

end


function ITensors.op!(Op::ITensor, ::OpName"RandU",
    ::SiteType"QCircuit", s::Index;
    α::Number, β::Number, γ::Number, δ::Number)

    Op[s'=>1, s=>1] = exp(1.0im*(α-β/2-δ/2))*cos(γ/2)
    Op[s'=>1, s=>2] = -exp(1.0im*(α-β/2+δ/2))*sin(γ/2)
    Op[s'=>2, s=>1] = exp(1.0im*(α+β/2-δ/2))*sin(γ/2)
    Op[s'=>2, s=>2] = exp(1.0im*(α+β/2+δ/2))*cos(γ/2)

    #U = exp(-1.0im * θ * [cos(α) sin(α)*cos(ϕ)-1.0im*sin(α)*sin(ϕ); 1.0im*sin(α)*sin(ϕ) -cos(α)])

end


################
# Two-site gates
################


function ITensors.op!(Op::ITensor, ::OpName"CNOT",
    ::SiteType"QCircuit", s::Index, t::Index)

    Op[s'=>1, s=>1, t'=>1, t =>1] = Complex(1.)
    Op[s'=>1, s=>1, t'=>2, t =>2] = Complex(1.)
    Op[s'=>2, s=>2, t'=>1, t =>2] = Complex(1.)
    Op[s'=>2, s=>2, t'=>2, t =>1] = Complex(1.)

end


function ITensors.op!(Op::ITensor, ::OpName"SWAP",
    ::SiteType"QCircuit", s::Index, t::Index)

    Op[s'=>1, s=>1, t'=>1, t =>1] = Complex(1.)
    Op[s'=>2, s=>1, t'=>1, t =>2] = Complex(1.)
    Op[s'=>1, s=>2, t'=>2, t =>1] = Complex(1.)
    Op[s'=>2, s=>2, t'=>2, t =>2] = Complex(1.)

end

"""
function ITensors.op!(Op::ITensor, ::OpName"CU",
    ::SiteType"QCircuit", s::Index; mat::Matrix{ComplexF64})

    Op[s'=>1, s=>1, t'=>1, t =>1] = Complex(1.)
    Op[s'=>1, s=>1, t'=>2, t =>2] = Complex(1.)
    Op[s'=>2, s=>2, t'=>1, t =>2] = mat[1,2]
    Op[s'=>2, s=>2, t'=>2, t =>1] = mat[2,1]
    Op[s'=>2, s=>2, t'=>1, t =>1] = mat[1,1]
    Op[s'=>2, s=>2, t'=>2, t =>2] = mat[2,2]
    #Op[s'=>1, s=>1] = U[1,1]
    #Op[s'=>1, s=>2] = U[1,2]
    #Op[s'=>2, s=>1] = U[2,1]
    #Op[s'=>2, s=>2] = U[2,2]

end
"""

function ITensors.op!(Op::ITensor, ::OpName"CU",
    ::SiteType"QCircuit", s::Index, t::Index;
    U11::Number, U12::Number, U21::Number, U22::Number)

    Op[s'=>1, s=>1, t'=>1, t =>1] = Complex(1.)
    Op[s'=>1, s=>1, t'=>2, t =>2] = Complex(1.)
    Op[s'=>2, s=>2, t'=>1, t =>2] = U11 #mat[1,2]
    Op[s'=>2, s=>2, t'=>2, t =>1] = U12 #mat[2,1]
    Op[s'=>2, s=>2, t'=>1, t =>1] = U21 #mat[1,1]
    Op[s'=>2, s=>2, t'=>2, t =>2] = U22 #mat[2,2]
    #Op[s'=>1, s=>1] = U[1,1]
    #Op[s'=>1, s=>2] = U[1,2]
    #Op[s'=>2, s=>1] = U[2,1]
    #Op[s'=>2, s=>2] = U[2,2]

end


# for quantum Fourier transform
function ITensors.op!(Op::ITensor, ::OpName"CRn",
    ::SiteType"QCircuit", s::Index, t::Index; n::Number)

    Op[s'=>1, s=>1, t'=>1, t =>1] = Complex(1.)
    Op[s'=>1, s=>1, t'=>2, t =>2] = Complex(1.)
    Op[s'=>2, s=>2, t'=>1, t =>2] = Complex(1.)
    #Op[s'=>2, s=>2, t'=>2, t =>1] = U12 #mat[2,1]
    #Op[s'=>2, s=>2, t'=>1, t =>1] = U21 #mat[1,1]
    Op[s'=>2, s=>2, t'=>2, t =>2] = exp(sign(n)*2π*1.0im/2^abs(n))

end


##################
# Three-site gates
##################


function ITensors.op!(Op::ITensor, ::OpName"TOFFOLI",
    ::SiteType"QCircuit", s::Index, t::Index, p::Index)

    Op[s'=>1, s=>1, t'=>1, t =>1, p'=>1, p=>1] = Complex(1.)
    Op[s'=>1, s=>1, t'=>1, t =>1, p'=>2, p=>2] = Complex(1.)
    Op[s'=>1, s=>1, t'=>2, t =>2, p'=>1, p=>1] = Complex(1.)
    Op[s'=>1, s=>1, t'=>2, t =>2, p'=>2, p=>2] = Complex(1.)
    Op[s'=>2, s=>2, t'=>1, t =>1, p'=>1, p=>1] = Complex(1.)
    Op[s'=>2, s=>2, t'=>1, t =>1, p'=>2, p=>2] = Complex(1.)
    Op[s'=>2, s=>2, t'=>2, t =>2, p'=>1, p=>2] = Complex(1.)
    Op[s'=>2, s=>2, t'=>2, t =>2, p'=>2, p=>1] = Complex(1.)

    #Op[s'=>2, s=>2, t'=>1, t =>2] = Complex(1.)
    #Op[s'=>2, s=>2, t'=>2, t =>1] = Complex(1.)

end



#############
"""
N = 4
θ = π
s = Index(2, "QCircuit")
X = op("X", s)
RX = op("Rz", s; θ)
Xmat = [0. 1.; 1. 0.]
println("X: ", array(X)==Xmat)
println("X: ", typeof(array(X)))


# get indices in QC Hilbert space, build MPO for single-site operator
sites = siteinds("QCircuit", N)
X_MPO = MPO(sites, ["E", "X", "E", "X"])
#RX_MPO = MPO(sites, ["E", "Rx", "E", "Rx"])
MPO_list = [op("E", sites[1]), op("Ry", sites[2]; θ), op("E", sites[3]), op("Ry", sites[4]; θ)]
Rz_MPO = MPO([op("E", sites[1]), op("Ry", sites[2]; θ), op("E", sites[3]), op("Ry", sites[4]; θ)])
Z_MPO = MPO(sites, ["E", "Z", "E", "Z"])
list1 = []
for i in 1:N
    push!(list1, "H")
end
H_MPO = MPO(sites, String.(list1))
println("X MPO: ", H_MPO)
println("Z MPO: ", H_MPO)
println("H MPO: ", H_MPO)
println("Rz MPO: ", Rz_MPO)

# add two MPOs and check that bond dimension has doubled
XpZ_MPO = add(X_MPO, Z_MPO)
println(XpZ_MPO)

# initialise state
psi = productMPS(sites, "0")
println("intial: ", psi)
println(ITensors.sample!(psi))

# apply MPO to state
psi = contract(Rz_MPO, psi)
println("after Rz: ", psi)
println(ITensors.sample!(psi))


os = [("CNOT", 4, 2)]
gate = ops(os, sites)
#@show os
println(gate)
psi = productMPS(sites, "1")
psi = apply(gate, psi)
println(ITensors.sample!(psi))

#################################

# initialise state
sites = siteinds("QCircuit", N)
psi = productMPS(sites, "0")
#println(psi)
println("initial: ", ITensors.sample!(psi))

# construct single site gate
X_gate = op("X", sites[1])
println(X_gate)

# apply single-site gate
newA = X_gate*psi[1]
noprime!(newA)
psi[1] = newA
println("after X: ", ITensors.sample!(psi))

# construct two-site gate
SWAP = op("SWAP", sites[1], sites[2])
#println(CNOT)

# apply two-site gate
orthogonalize!(psi, 1)
wf = (psi[1]*psi[2])*SWAP
noprime!(wf)
ind = uniqueinds(psi[1], psi[2])
println("check indices: ", ind)
U, S, V = svd(wf, ind, cutoff=1E-8)
println("U ", U)
println("S ", S)
println("V ", V)
psi[1] = U
psi[2] = S*V
println(psi)

#println(maxlinkdim(psi))
s = ITensors.sample(orthogonalize!(psi, 1))
#s = ITensors.sample(psi)
println("final: ", s)


"""
