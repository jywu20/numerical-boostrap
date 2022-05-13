#region how verbose the program is 

#const working_path = "./"
#const output_name = "2022-4-30-run-1-res"

const working_path = "D:\\Projects\\numerical-boostrap\\hubbard-2d-prototype\\"
const output_name = "2022-5-13-run-1-res.out"

const full_output_name = working_path * output_name

# If there exists working_path * output_name already, throw an error
const no_conflict = false

# Display operator labels 
const show_operator_labels = true

# Display operators involved in the bootstrap process
const show_hubbard_opstr_basis = true 

# Display the M matrix 
const show_M_matrix = true

# Display the constraints
const show_constraints = true

# Display the Hamiltonian
const show_hamiltonian = true

#endregion

#region model parameter and cutoff

# The hopping parameter and local repulsion  
U = 4.0
t = 1.0

# l(O) ≤ K cutoff
K = 10
site_num = (2K + 1)^2

# When this flag is `true`, no actual optimization will be done. For debugging only.
no_optimization = false 
# Feasibility report 
feasibility_check = false 

no_geometric_symmetry = false

max_iter = 10000

#endregion