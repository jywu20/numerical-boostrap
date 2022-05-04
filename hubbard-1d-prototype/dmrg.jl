using ITensors

E, ψ = let 
    N = 100
    bond_num = 10

    sites = siteinds("Electron", N; conserve_qns = true)

    ham_ampo = OpSum()

    for i = 1 : N
        ham_ampo += U ,"Nupdn", i
    end

    for i in 1 : N - 1
        ham_ampo += - t, "Cdagup", i, "Cup", i + 1
        ham_ampo += - t, "Cdagdn", i, "Cdn", i + 1
        ham_ampo += - t, "Cdagup", i + 1, "Cup", i
        ham_ampo += - t, "Cdagdn", i + 1, "Cdn", i
    end

    H = MPO(ham_ampo, sites)

    state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    ψ0 = productMPS(sites,state)

    @show flux(ψ0)

    sweeps = Sweeps(5)
    setmaxdim!(sweeps, 10,20,100,100,200)
    setcutoff!(sweeps, 1E-10)

    dmrg(H, ψ0, sweeps)
end

