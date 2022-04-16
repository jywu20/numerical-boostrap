using JuMP, COSMO 

model = Model(COSMO.Optimizer)
set_optimizer_attributes(model, "max_iter" => max_iter)

#region Create variables corresponding to the real and imaginary parts of correlation functions
# We put the real and imaginary part of `hubbard_opstr_basis` apart, because COSMO only supports 
# real variables
@variable(model, hubbard_opstr_basis_expected[1 : 2hubbard_opstr_basis_length])

"""
- `i`: the indices in `hubbard_opstr_basis`
"""
function hubbard_opstr_basis_expected_real_imag_parts(i, real_or_imag) 
    if real_or_imag == :real
        return hubbard_opstr_basis_expected[2i - 1]
    end
    if real_or_imag == :imag
        return hubbard_opstr_basis_expected[2i]
    end
end

I22 = Matrix(1*I, 2, 2)
O22 = 0.0 * I22
Im22 = [
    0   -1 ; 
    1    0
]

hubbard_opstr_basis_expected_real = [hubbard_opstr_basis_expected_real_imag_parts(i, :real) * I22 
    for i in 1 : hubbard_opstr_basis_length]
hubbard_opstr_basis_expected_imag = [hubbard_opstr_basis_expected_real_imag_parts(i, :imag) * Im22 
    for i in 1 : hubbard_opstr_basis_length]

function complex_to_mat(coefficients)
    real_part = transpose(real(coefficients))
    imag_part = transpose(imag(coefficients))
    real_part_mat_version = map(x -> x * I22, real_part)
    imag_part_mat_version = map(x -> x * Im22, imag_part)
    (real_part_mat_version + imag_part_mat_version) * (hubbard_opstr_basis_expected_real + hubbard_opstr_basis_expected_imag)
end

#endregion

#region Imposing constraints

for constraint_coefficients in H_constraints_coefficients
    @constraint(model, complex_to_mat(constraint_coefficients) .== O22)
end

@variable(model, 
    M[1 : 2 * length(M_mat_spanning_opstr_indices), 1 : 2 * length(M_mat_spanning_opstr_indices)], PSD)

# Here i and j are indices of O_i and O_j which define M_{ij} = ⟨O_i O_j⟩
for i in 1 : length(M_mat_spanning_opstr_indices)
    for j in 1 : length(M_mat_spanning_opstr_indices)
        @constraint(model, M[2i - 1 : 2i, 2j - 1 : 2j] .== complex_to_mat(M_coefficient[i, j]))
    end
end

# Use the Hamiltonian on a single site as the optimization target
H_1 = - t * sum(map(Iterators.product([2, 4, 6, 8], [↑, ↓])) do label
    i = label[1]
    σ = label[2]
    cdag(1, σ) * c(i, σ)
end) + U * cdag(1, ↑) * c(1, ↑) * cdag(1, ↓) * c(1, ↓)

H_1 = normal_form(H_1)

# hubbard_opstr_coefficients(H_1) is supposed to be real, so we take its real part
@objective(model, Min, complex_to_mat(hubbard_opstr_coefficients(H_1))[1, 1])

#endregion