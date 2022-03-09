using JuMP, CSDP
using LinearAlgebra  
using OffsetArrays
using Base.Iterators
using Test

# Interaction strength
g = 1.0
# Maximal length of x operator sequence and p operator sequence ever considered in this program.
# It should be greater than 4 + the length of x/p sequence of O in ⟨[H, O]⟩=0, 
# so that when computing commutation relation with the Hamiltonian,
# there will be no out of bound error. Similarly, when constructing the M matrix, we need to 
# make sure that 2K ≤ L.
L_max = 12
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

@testset "Operator string labeling" begin
    @test [xpopstr_index(i, j) for i in xpoppower_range for j in xpoppower_range] == 0 : xpopspace_dim 
end

##

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

@testset "Operator string generating and stringifying." begin
    op0 = xpopstr_const(3.4)
    @test xpopstr_stringify(op0) == "3.4"
    op1 = xpopstr_xp_power(2, 3)
    @test xpopstr_stringify(op1) == "x^2 p^3"
    @test xpopstr_stringify(- op1) == "- x^2 p^3"
    @test xpopstr_stringify(2.0 * op1) == "2.0 x^2 p^3"
    @test xpopstr_stringify(- op0 + op1) == "- 3.4 + x^2 p^3"
    op2 = xpopstr_xp_power(4, 3)
    @test xpopstr_stringify(- op0 + 2.3 * op1 - im * op2) == "- 3.4 + 2.3 x^2 p^3 - 1.0im x^4 p^3"
    @test xpopstr_stringify(- op0 + 2.3 * op1 + (1 + im) * op2) == "- 3.4 + 2.3 x^2 p^3 + (1.0 + 1.0im) x^4 p^3"
    @test xpopstr_stringify(- op0 + 2.3 * op1 + (- 2.3 + im) * op2) == "- 3.4 + 2.3 x^2 p^3 + (-2.3 + 1.0im) x^4 p^3"
    op3 = xpopstr_xp_power(0, 3)
    @test xpopstr_stringify(im * op3) == "1.0im p^3"
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

@testset "Operator multiplication by index offset" begin
    @test xpopstr_left_x_right_p_mul(2, 2, 1) == 2 + xpopstr_left_x_right_p_mul_offset(2, 1)
    @test xpopstr_left_x_right_p_mul(10, 2, 3) == 10 + xpopstr_left_x_right_p_mul_offset(2, 3)
end

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

@testset "Multiplication in the form of x^x_power * O * p^p_power" begin
    op0 = xpopstr_const(3.4)
    @test xpopstr_stringify(xpopstr_left_x_right_p_mul(op0, 2, 3)) == "3.4 x^2 p^3"
    op1 = xpopstr_const(2.5) + im * xpopstr_xp_power(3, 4)
    @test xpopstr_stringify(xpopstr_left_x_right_p_mul(op1, 4, 2)) == "2.5 x^4 p^2 + 1.0im x^7 p^6"
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

xpopstr_comm(x_power_1, p_power_1, x_power_2, p_power_2) = xpopstr_normal_ord(x_power_1, p_power_1, x_power_2, p_power_2) - xpopstr_normal_ord(x_power_2, p_power_2, x_power_1, p_power_1)

@testset "Normal ordering of x^x_power_1 p^p_power_1 x^x_power_2 p^p_power_2" begin
    # See ./commutation-x-p.nb for symbolic benchmarks 
    @test xpopstr_stringify(xpopstr_normal_ord(2, 3, 3, 2)) == 
        "6.0im x^2 p^2 - 18.0 x^3 p^3 - 9.0im x^4 p^4 + x^5 p^5"
    @test xpopstr_stringify(xpopstr_normal_ord(2, 3, 1, 3)) == "- 3.0im x^2 p^5 + x^3 p^6"
end

@testset "Commutation between x^x_power_1 p^p_power_1 and x^x_power_2 p^p_power_2" begin
    # See ./commutation-x-p.nb for symbolic benchmarks 
    @test xpopstr_stringify(xpopstr_comm(2, 3, 1, 3)) == "6.0 x p^4 + 3.0im x^2 p^5" 
    @test xpopstr_stringify(xpopstr_comm(2, 3, 2, 2)) == "- 4.0 x^2 p^3 - 2.0im x^3 p^4" 
end

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

@testset "Commutator with the Hamiltonian" begin
    op1 = xpopstr_xp_power(2, 1)
    @test xpopstr_stringify(comm_with_ham(op1)) == "2.0 p + 4.0im x p^2 - 2.0im x^3 - 4.0im x^5"
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

@testset "Finding the maximal power of x and p in an operator" begin
    op0 = xpopstr_const(3.4)
    op1 = xpopstr_const(2.5) + im * xpopstr_xp_power(3, 4)
    @test max_x_power(op0) == 0
    @test max_p_power(op0) == 0
    @test max_x_power(op1) == 3
    @test max_p_power(op1) == 4
end

##

model = Model(CSDP.Optimizer)
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

variable_list_real = [xpopstr_expected_real_imag_parts(i, :real) * I22 for i in 1 : xpopspace_dim]
variable_list_imag = [xpopstr_expected_real_imag_parts(i, :imag) * Im22 for i in 1 : xpopspace_dim]
xpopstr_basis_real = OffsetArray([I22, variable_list_real...], xpopspace_index_range)
xpopstr_basis_imag = OffsetArray([O22, variable_list_imag...], xpopspace_index_range)

# Building constraints
# make sure the size of OH does not cause any out of bound error.
# Then imposing constraints to the variables

for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    lhs = transpose(real(cons)) * xpopstr_basis_real + transpose(imag(cons)) * xpopstr_basis_imag
    @constraint(model, lhs .== O22)
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
        op_ij = xpopstr_normal_ord(op1_idx_xpower, op1_idx_ppower, op2_idx_xpower, op2_idx_ppower)

        real_part = transpose(real(op_ij)) * xpopstr_basis_real * I22
        imag_part = transpose(imag(op_ij)) * xpopstr_basis_imag * Im22
        @constraint(model, M[2i - 1 : 2i, 2j - 1 : 2j] .== real_part + imag_part)
    end
end

@objective(model, Min, xpopstr_expected[xpopstr_index(2, 0)] + xpopstr_expected[xpopstr_index(0, 2)] + g * xpopstr_expected[xpopstr_index(4, 0)])
optimize!(model)

@show objective_value(model)
