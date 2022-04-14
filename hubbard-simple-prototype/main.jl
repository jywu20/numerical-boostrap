#region Start up
# Load configuration
include("configuration.jl")

# If there exists full_output_name already, hint the user
if isfile(full_output_name)
    if no_conflict
        throw("Error: output file $(full_output_name) already exists.")
    end

    @warn "Output file $(full_output_name) already exists. Its content will be overwritten."
end

open(full_output_name, "w") do file
    println(file, "2D square lattice Hubbard model bootstrap")
    println(file, "Configuration loaded.")
    println(file, "")
end
#endregion

#region Generate labels for operators (i.e. "site 34, spin ↑")
include("operator-label.jl")

open(full_output_name, "a") do file
    println(file, "Operator labels constructed.")
end
#endregion

#region Find M matrix and 
include("operator-algebra.jl")

open(full_output_name, "a") do file
    println(file, "Constraints and M matrix assembled.")
end
#endregion

#region Construct the Optimization problem

include("optimization_problem.jl")

open(full_output_name, "a") do file
    println(file, "Optimization problem constructed. Starting optimizing.")
    println(file, "--------------------------------------")
    println(file)
end

include("run_optimization.jl")

include("correlation-functions.jl")

#endregion
