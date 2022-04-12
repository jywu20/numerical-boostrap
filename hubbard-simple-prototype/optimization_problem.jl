using JuMP, COSMO 

model = Model(COSMO.Optimizer)
set_optimizer_attributes(model, "max_iter" => max_iter)

# We put the real and imaginary part of `hubbard_opstr_basis` apart, because COSMO only supports 
# real variables
@variable(model, hubbard_opstr_basis_expected[1 : 2hubbard_opstr_basis_length])

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

for constraint_coefficients in H_constraints_coefficients
    
end