if ! no_optimization
    optimize!(model)

    open(full_output_name, "a") do file
        println(file, "-----------------------------------------------------------")
        println(file, "Optimization results:")
        println(file, "")
        println(file, "Objective value:   $(objective_value(model))")
    
        println(file)
        println(file, "-----------------------------------------------------------")
        println(file)
    end
end
