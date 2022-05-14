final_expected_real_value(op) = op |> normal_form |> hubbard_opstr_coefficients |> coefficients_to_variable_ref

if ! no_optimization
    solution_point = map(all_variables(model)) do optimize_var
        optimize_var => value(optimize_var)
    end |> Dict

    if dmrg_benchmark
        open(full_output_name, "a") do file
            println(file, "variable  DMRG benchmark value    bootstrap value    relative error ")
            for optimize_var in keys(jump_benchmark_point)
                dmrg_value = jump_benchmark_point[optimize_var]
                bootstrap_value = solution_point[optimize_var]
                relative_err = abs((dmrg_value - bootstrap_value) / max(abs(dmrg_value), abs(bootstrap_value)))
                println(file, "$optimize_var    $dmrg_value    $bootstrap_value     $relative_err")
            end
            println(file)
            println(file)
        end
    end

    open(full_output_name, "a") do file
        #region ⟨n_i⟩
    
        n_i_var = (cdag(1, ↑) * c(1, ↑) * cdag(1, ↓) * c(1, ↓)) |> final_expected_real_value
        println(file, "Expected ⟨n_i↑ n_i↓⟩:  $(value(n_i_var))")
        println(file)
    
        #endregion
    
        #region ⟨S_{zi}⟩
        Sz_i_var = (cdag(1, ↑) * c(1, ↑) - cdag(1, ↓) * c(1, ↓)) |> final_expected_real_value
        println(file, "Expected ⟨Sz_i⟩:  $(value(Sz_i_var))")
        println(file)
        #endregion
    end
end