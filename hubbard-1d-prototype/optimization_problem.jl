using JuMP, COSMO 

model = Model(COSMO.Optimizer)
set_optimizer_attributes(model, "max_iter" => max_iter)

#region Create variables corresponding to the real and imaginary parts of correlation functions
# We put the real and imaginary part of `hubbard_opstr_basis` apart, because COSMO only supports 
# real variables

hubbard_opstr_basis_expected = Vector(undef, hubbard_opstr_basis_length)
hubbard_opstr_basis_expected[1] = 1.0
for opstr_idx in 2 : hubbard_opstr_basis_length
    op = hubbard_opstr_basis[opstr_idx]
    op_string_form = string(op)
    if ! (op in particle_number_constraint_ops) && ! (op in spin_constraint_ops)
        hubbard_opstr_basis_expected[opstr_idx] = @variable(model, base_name = "⟨$op_string_form⟩")
    else
        hubbard_opstr_basis_expected[opstr_idx] = 0.0
    end
end

function coefficients_to_variable_ref(constraint::AbstractVector)
    res = 0
    for (i, c) in enumerate(constraint)
        res += c * hubbard_opstr_basis_expected[i]
    end
    res
end

if no_geometric_symmetry
    all_constraints = H_constraints_coefficients
else
    all_constraints = [
        H_constraints_coefficients..., 
        translational_constraint_coefficients..., 
        reflectional_constraints_coefficients...]
end

for constraint_coefficients in all_constraints
    if constraint_coefficients != hubbard_opstr_zero
        @constraint(model, coefficients_to_variable_ref(constraint_coefficients) == 0)
    end
end

for i in 1 : length(site_list)
    if ! (cdag(i, ↑) * c(i, ↑) in hubbard_opstr_basis)
        break
    end
    total_particle_number_constraint = coefficients_to_variable_ref(hubbard_opstr_coefficients(cdag(i, ↑) * c(i, ↑) + cdag(i, ↓) * c(i, ↓)))
    @constraint(model, total_particle_number_constraint == 1)
end

if ! feasibility_check 
    @variable(model, 
        M[1 : length(M_mat_spanning_opstr_indices), 1 : length(M_mat_spanning_opstr_indices)], PSD)
else
    @variable(model, 
        M[1 : length(M_mat_spanning_opstr_indices), 1 : length(M_mat_spanning_opstr_indices)])
end

# Here i and j are indices of O_i and O_j which define M_{ij} = ⟨O_i O_j⟩
for i in 1 : length(M_mat_spanning_opstr_indices)
    for j in 1 : length(M_mat_spanning_opstr_indices)
        @constraint(model, M[i, j] == coefficients_to_variable_ref(M_coefficient[i, j]))
    end
end

# Use the Hamiltonian on a single site as the optimization target
H_1 = - t * sum(map(Iterators.product([2, 3], [↑, ↓])) do label
    i = label[1]
    σ = label[2]
    cdag(1, σ) * c(i, σ)
end) + U * cdag(1, ↑) * c(1, ↑) * cdag(1, ↓) * c(1, ↓)

H_1 = normal_form(H_1)

# hubbard_opstr_coefficients(H_1) is supposed to be real, so we take its real part
@objective(model, Min, coefficients_to_variable_ref(hubbard_opstr_coefficients(H_1)))

#endregion