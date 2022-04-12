#region how verbose the program is 

const working_path = "D:\\Projects\\numerical-boostrap\\hubbard-simple-prototype\\"
const output_name = "result"

const full_output_name = working_path * output_name

# If there exists working_path * output_name already, throw an error
const no_conflict = false

# Display operators involved in the bootstrap process
const show_hubbard_opstr_basis = true

#endregion

#region model parameter and cutoff

# The hopping parameter and local repulsion  
U = 4.0
t = 1.0

# l(O) ≤ K cutoff
K = 5
site_num = (2K + 1)^2

max_iter = 10000

#endregion