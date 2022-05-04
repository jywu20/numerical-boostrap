using ITensors

N = 100

E, ψ, sites = let 
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
    setmaxdim!(sweeps, 10,20,100,200,200)
    setcutoff!(sweeps, 1e-10)

    E, ψ = dmrg(H, ψ0, sweeps)
    E, ψ, sites
end

##

benchmark_point_dmrg = Dict{QuExpr, Float64}()

benchmark_point_dmrg[hubbard_opstr_basis[1]] = 1.0

let 
    sites_middle = Int(N / 2)

    for op in hubbard_opstr_basis[2 : end]
        opstr_string = string(op)
        opstr_string = replace(opstr_string, "c" => "C")
        opstr_string = replace(opstr_string, "†" => "dag")
        opstr_each_site_strings = split(opstr_string)

        Sz = 0
        Nf = 0

        opstr_each_site_op_and_idx = map(opstr_each_site_strings) do str
            str = str[1 : end - 2]

            if str[end] == '-'
                str = str[1 : end - 1]
                head, site_idx = split(str, "(")
                str = head * "dn"
                Szi = -1
            else 
                head, site_idx = split(str, "(")
                str = head * "up"
                Szi = 1
            end

            if head == "Cdag"
                Nf += 1
                Sz += Szi
            else 
                Nf -= 1
                Sz -= Szi
            end

            [str, sites_middle + site_list[parse(Int, site_idx)]]
        end

        if isodd(length(opstr_each_site_op_and_idx))
            push!(benchmark_point_dmrg, op => 0.0)
        elseif Sz == 0 && Nf == 0
            opstr_def = vcat(opstr_each_site_op_and_idx...) |> Tuple

            opstr_mpo = OpSum()
            opstr_mpo += opstr_def

            push!(benchmark_point_dmrg, op => inner(ψ, MPO(opstr_mpo, sites), ψ))
        else 
            push!(benchmark_point_dmrg, op => 0.0)
        end
    end
end