using ITensors

function test_MPS()
  N = 3
  chi = 4
  sites = siteinds("S=1/2",N)

  A1=[1 0]
  A2=[0 1]
  A3=[1/2 1/2]
  A=kron(A1,A2,A3)
  A=A/norm(A)
  print(A,"\n")
  psi = MPS(A,sites;cutoff=1e-8,maxdim=10)
  print(psi[1],"\n")
  print(psi[2],"\n")
  println(psi[3])
  print(psi[1]*psi[2]*psi[3],"\n")
  magz = ITensors.expect(psi,"Sz")
  for (j,mz) in enumerate(magz)
    println("$j $mz")
  end
end

function test_MPO()
  N = 10
  sites = siteinds("S=1/2",N)
  #sites2 = siteinds("S=1/2",N)

  ampo = OpSum()
  for j=1:N-1
    ampo += "Sz",j,"Sz",j+1
    ampo += 1/2,"S+",j,"S-",j+1
    ampo += 1/2,"S-",j,"S+",j+1
  end
  H = MPO(ampo,sites)
  psi0 = randomMPS(sites,10)
  psi2 = randomMPS(sites,20)
  psi1 = contract(H,psi0,maxdim=10)
  println("works")
  println(inner(psi1,psi0))
  println("works2")

  println(psi0)
  #println(psi1)
  println(psi1)
  #deepcopy(psi0) + deepcopy(psi1)
  add(psi0, psi2)
  add(psi0, psi1)
  #println(psi0+psi1)
  #println(add(psi0,psi1))
end

test_MPO()

#psi0 = randomMPS(sites,10)
#psi2 = randomMPS(sites,20)
#add(psi0, psi2)


#O1=[1 0; 0 -1]
#id=[1 0; 0 1]
#h=kron(id,id,id)
#println(h*(transpose(A)))
#H=MPO(h,sites;cutoff=1e-8,maxdim=10)

#psi_2 = contract(psi,H)
#f = inner(psi_2,psi)
#print(f,"\n")
#print(inner(psi,psi),"\n")
#for (j,mz) in enumerate(magz)
#    println("$j $mz")
#end
