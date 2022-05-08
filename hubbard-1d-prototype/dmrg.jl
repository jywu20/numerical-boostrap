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

#region 

jump_benchmark_point = Dict{VariableRef, Float64}()
for (opstr_idx, op) in enumerate(hubbard_opstr_basis)
    optimization_var = hubbard_opstr_basis_expected[opstr_idx]
    if typeof(optimization_var) == Float64
        # println(optimization_var, " ", benchmark_point_dmrg[op])
        continue
    end
    push!(jump_benchmark_point, optimization_var => benchmark_point_dmrg[op])
end

for (i, opstr_index_1) in enumerate(M_mat_spanning_opstr_indices)
    for (j, opstr_index_2) in enumerate(M_mat_spanning_opstr_indices)
        opstr_1 = hubbard_opstr_basis[opstr_index_1]
        opstr_2 = hubbard_opstr_basis[opstr_index_2]

        Mij = normal_form(opstr_1' * opstr_2)
        
        res = 0.0
        for (term_label, coefficient) in Mij.terms
            term_label = QuExpr(Dict(term_label => 1.0))
            res += benchmark_point_dmrg[term_label] * coefficient
        end

        push!(jump_benchmark_point, M[i, j] => res)
    end
end

if feasibility_check
    feasibility_report = primal_feasibility_report(model, jump_benchmark_point)

    M_dmrg = map(M) do optimize_var
        jump_benchmark_point[optimize_var]
    end

    M_dmrg_spectrum = eigen(M_dmrg).values
end

function dmrg_value(op)
    value(op |> hubbard_opstr_coefficients |> coefficients_to_variable_ref) do optimize_var
        jump_benchmark_point[optimize_var]
    end
end

#endregion