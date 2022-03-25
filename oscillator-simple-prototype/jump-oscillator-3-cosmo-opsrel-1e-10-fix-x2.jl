using JuMP, COSMO 
using LinearAlgebra  
using OffsetArrays
using Base.Iterators

# Interaction strength
g = 1.0
# Maximal length of x operator sequence and p operator sequence ever considered in this program.
# It should be greater than 4 + the length of x/p sequence of O in ⟨[H, O]⟩=0, 
# so that when computing commutation relation with the Hamiltonian,
# there will be no out of bound error. Similarly, when constructing the M matrix, we need to 
# make sure that 2K ≤ L.
L_max = 5 
# The dimension of the operator space; the -1 term comes from the fact that a constant is not 
# considered as an operator 
xpopspace_dim = (2L_max + 1)^2 - 1
# The range of possible x / p power
xpoppower_range = 0 : 2L_max
# The complete range of indexes of operators, the length of which is xpopspace_dim + 1 = (2L_max + 1)^2
xpopspace_index_range = 0 : xpopspace_dim

# Define variables. We label the expectation values as the follows: x^m p^n is labeled as m * (K + 1) + n 
# Note: the index starts from 0 and should be used together with an OffsetArray
xpopstr_index(x_power, p_power) = x_power * (2L_max + 1) + p_power
index_to_xpower(idx) = Int(floor(idx / (2L_max + 1)))
index_to_ppower(idx) = idx % (2L_max + 1)

"""
The coefficient list of a constant "operator" c.
"""
function xpopstr_const(c)
    result = OffsetArray(zeros(ComplexF64, xpopspace_dim + 1), xpopspace_index_range)
    result[0] = c
    result
end

"""
x^x_power p^p_power
"""
function xpopstr_xp_power(x_power, p_power)
    result = xpopstr_const(0.0)
    result[xpopstr_index(x_power, p_power)] = 1.0
    result
end

function power_str(var, pow)
    if pow == 0
        return ""
    end
    if pow == 1
        return var 
    end
    return "$var^$pow"
end

function xpopstr_power_str(x_power, p_power)
    if x_power != 0 && p_power != 0
        return " $(power_str("x", x_power)) $(power_str("p", p_power))"
    end
    if x_power == 0 && p_power != 0
        return " $(power_str("p", p_power))"
    end
    if x_power != 0 && p_power == 0
        return " $(power_str("x", x_power))"
    end
    return ""
end

function coefficient_format(c)
    if c == 1.0
        return " +"
    end
    if c == - 1.0
        return " -"
    end
    if imag(c) == 0
        c = real(c)
        if c > 0 
            return " + $c"
        end
        if c < 0
            return " - $(abs(c))"
        end
    end
    if real(c) == 0
        c = imag(c)
        if c > 0
            return " + $(c)im"
        end
        if c < 0
            return " - $(- c)im"
        end
    end
    " + ($(string(c)))"
end

function xpopstr_stringify(coefficients::OffsetArray)
    result_terms = []
    for idx in xpopspace_index_range
        x_power = index_to_xpower(idx)
        p_power = index_to_ppower(idx)
        if coefficients[idx] != 0
            push!(result_terms, 
                "$(coefficient_format(coefficients[idx]))$(xpopstr_power_str(x_power, p_power))")
        end
    end
    output_str = strip(join(result_terms, ""))
    if output_str[1] == '+'
        return strip(output_str[2:end])
    end
    return output_str
end

"""
x^x_power * idx * p^p_power
"""
function xpopstr_left_x_right_p_mul(idx, x_power, p_power)
    x_power_0 = index_to_xpower(idx)
    p_power_0 = index_to_ppower(idx)
    x_power′ = x_power_0 + x_power
    p_power′ = p_power_0 + p_power
    xpopstr_index(x_power′, p_power′)
end

xpopstr_left_x_right_p_mul_offset(x_power, p_power) = x_power * (2L_max + 1) + p_power

"""
x^x_power * coefficients * p^p_power
May throw out of index error if the highest order term is beyond xpopspace_index_range
"""
function xpopstr_left_x_right_p_mul(coefficients, x_power, p_power)
    offset = xpopstr_left_x_right_p_mul_offset(x_power, p_power)
    result = xpopstr_const(0.0)
    for idx in offset : xpopspace_dim
        result[idx] = coefficients[idx - offset]
    end
    result
end

"""
Calculate coefficients of [x^n, p^m] according to the McCoy's formula (not the operator, 
because this function is used to calculate the normal ordering)
"""
function xpopstr_comm(x_power, p_power)
    n = x_power
    m = p_power
    result = xpopstr_const(0.0)
    for k in 1 : min(m, n)
        numerator = - (- im)^k * factorial(n) * factorial(m)
        denominator = factorial(k) * factorial(n - k) * factorial(m - k)
        result[xpopstr_index(n - k, m - k)] = numerator / denominator
    end
    result
end

"""
Normal ordered version of x^x_power_1 p^p_power_1 x^x_power_2 p^p_power_2
"""
function xpopstr_normal_ord(x_power_1, p_power_1, x_power_2, p_power_2)
    # We name the first term on the RHS of
    # x^x_power_1 p^p_power_1 x^x_power_2 p^p_power_2 = 
    # x^x_power_1 x^x_power_2 p^p_power_1 p^p_power_2 - x^x_power_1 [x^x_power_2, p^p_power_1] p^p_power_2
    # as pre_normal_ord_term, the second term as comm_term
    pre_normal_ord_term = xpopstr_xp_power(x_power_1 + x_power_2, p_power_1 + p_power_2)
    comm_term = xpopstr_left_x_right_p_mul(xpopstr_comm(x_power_2, p_power_1), x_power_1, p_power_2)
    pre_normal_ord_term - comm_term
end

xpopstr_comm(x_power_1, p_power_1, x_power_2, p_power_2) = 
    xpopstr_normal_ord(x_power_1, p_power_1, x_power_2, p_power_2) - xpopstr_normal_ord(x_power_2, p_power_2, x_power_1, p_power_1)

"""
For simplicity, we do not implement a full version of commutation; we just calculate the commutator between 
each term of an operator and the three terms of the Hamiltonian.
"""
function comm_with_ham(op::OffsetArray)
    result = xpopstr_const(0.0)
    for idx in xpopspace_index_range
        if op[idx] != 0
            x_power_in_op = index_to_xpower(idx)
            p_power_in_op = index_to_ppower(idx)

            # [op, x^2]
            result += xpopstr_comm(x_power_in_op, p_power_in_op, 2, 0)
            # [op, p^2]
            result += xpopstr_comm(x_power_in_op, p_power_in_op, 0, 2)
            # [op, g x^4]
            result += g * xpopstr_comm(x_power_in_op, p_power_in_op, 4, 0)
        end
    end
    result
end

function max_x_power(op::OffsetArray)
    nonzero_terms_idx = xpopspace_index_range[collect(op .!= 0)]
    if length(nonzero_terms_idx) == 0
        return 0
    end
    max(map(index_to_xpower, nonzero_terms_idx)...)
end

function max_p_power(op::OffsetArray)
    nonzero_terms_idx = xpopspace_index_range[collect(op .!= 0)]
    if length(nonzero_terms_idx) == 0
        return 0
    end
    max(map(index_to_ppower, nonzero_terms_idx)...)
end

model = Model(COSMO.Optimizer)
set_optimizer_attributes(model, "max_iter" => 30000000, "eps_rel" => 1.0e-10)
# Only non-constant operators have uncertain expectations
# Note: since a generic O is not Hermitian, we need to replace O, O† by (O + O†), i (O - O†)
# Note that we need to record both the imaginary part and the real part of each ⟨O⟩, so the 
# number of variables is doubled
@variable(model, xpopstr_expected[1 : 2xpopspace_dim])

xpopstr_expected_real_imag_parts(i, real_or_imag) = begin
    if real_or_imag == :real
        return xpopstr_expected[2i - 1]
    end
    if real_or_imag == :imag
        return xpopstr_expected[2i]
    end
end

I22 = Matrix(1*I, 2, 2)
O22 = 0.0 * I22
Im22 = [
    0   -1 ; 
    1    0
]

zero_xpopstr = xpopstr_const(0.0)

variable_list_real = [xpopstr_expected_real_imag_parts(i, :real) * I22 for i in 1 : xpopspace_dim]
variable_list_imag = [xpopstr_expected_real_imag_parts(i, :imag) * Im22 for i in 1 : xpopspace_dim]
xpopstr_basis_real = OffsetArray([I22, variable_list_real...], xpopspace_index_range)
xpopstr_basis_imag = OffsetArray([O22, variable_list_imag...], xpopspace_index_range)

function complex_to_mat(coefficients)
    real_part = transpose(real(coefficients))
    imag_part = transpose(imag(coefficients))
    real_part_mat_version = map(x -> x * I22, real_part)
    imag_part_mat_version = map(x -> x * Im22, imag_part)
    (real_part_mat_version + imag_part_mat_version) * (xpopstr_basis_real + xpopstr_basis_imag)
end

# Building constraints
# make sure the size of OH does not cause any out of bound error.
# Then imposing constraints to the variables

for x_power in 0 : 2L_max - 4, p_power in 0 : 2L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = real(cons)
    cons_imag = imag(cons)
    lhs = complex_to_mat(cons)
    if cons_real == cons_imag == zero_xpopstr
        continue
    end
    if cons_real == zero_xpopstr
        @constraint(model, lhs[1, 2] == 0.0)
    elseif cons_imag == zero_xpopstr
        @constraint(model, lhs[1, 1] == 0.0)
    else
        @constraint(model, lhs .== O22)
    end
end

# Construct the M matrix and impose the semidefinite constraint 
# The range of two indices of M matrix defined in papers is from 0 to L; note that this is NOT the index in the 
# JuMP constraint @variable(model, M[...]), because 1. we need to replace 0 with I and i with [0 -1; 1 0], 
# so the size of the PSD matrix is twice as much as the size of the M matrix in analytical calculation, and
# 2. in JuMP we count from 1 not 0 and OffsetArray cannot be created using @variable

M_index_to_xpopstr_index = zeros(Int, (L_max + 1)^2)
for i in 0 : L_max
    for j in 0 : L_max
        M_idx = i * (L_max + 1) + j + 1
        M_index_to_xpopstr_index[M_idx] = xpopstr_index(i, j)
    end
end

# To check if the opeartors are correct, run
# map(idx -> xpopstr_power_str(index_to_xpower(idx), index_to_ppower(idx)), M_index_to_xpopstr_index) 

@variable(model, M[1 : 2 * (L_max + 1)^2, 1 : 2 * (L_max + 1)^2], PSD)

for i in 1 : (L_max + 1)^2
    for j in i : (L_max + 1)^2
        op1_idx = M_index_to_xpopstr_index[i]
        op2_idx = M_index_to_xpopstr_index[j]
        op1_idx_xpower = index_to_xpower(op1_idx)
        op1_idx_ppower = index_to_ppower(op1_idx)
        op2_idx_xpower = index_to_xpower(op2_idx)
        op2_idx_ppower = index_to_ppower(op2_idx)
        op_ij = xpopstr_normal_ord(0, op1_idx_ppower, op1_idx_xpower + op2_idx_xpower, op2_idx_ppower)

        @constraint(model, M[2i - 1 : 2i, 2j - 1 : 2j] .== complex_to_mat(op_ij))
    end
end

@objective(model, Min, 
    xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :real) + 
    xpopstr_expected_real_imag_parts(xpopstr_index(0, 2), :real) + 
    g * xpopstr_expected_real_imag_parts(xpopstr_index(4, 0), :real))

# The initial value of E and ⟨x²⟩ are [1.37, 0.298], and therefore 
# x⁴ = (E - 2x²) / 3g = 0.258 
# x²_0 = 0.298
# x⁴_0 = 0.258
# set_start_value(xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :real), x²_0)
# set_start_value(xpopstr_expected_real_imag_parts(xpopstr_index(4, 0), :real), x⁴_0)

# Fix one varialbe and see what's going on
# 
x²_0_standard = 0.298
fix(xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :real), x²_0_standard)
fix(xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :imag), 0.0)

optimize!(model)

println("-----------------------------------------------------------")
println("Results:")
println("")
println("Objective value:   $(objective_value(model))")
println("x square expectation:     $(value(xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :real)))")
