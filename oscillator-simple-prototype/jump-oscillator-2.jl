using JuMP, CSDP
using LinearAlgebra  
using OffsetArrays
using Test

# Interaction strength
g = 1.0
# Maximal length of x operator sequence and p operator sequence in C_2
L_max = 7
# The dimension of the operator space; the -1 term comes from the fact that a constant is not 
# considered as an operator 
xpopspace_dim = (2L_max + 1)^2 - 1
# The range of possible x / p power 
xpoppower_range = 0 : 2L_max
# The complete range of indexes of operators 
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

model = Model(CSDP.Optimizer)
# Only non-constant operators have uncertain expectations
@variable(model, xpopstr_expected[1 : xpopspace_dim])

"""
The coefficient list of a constant "operator" c.
"""
function xpopstr_const(c)
    result = OffsetArray(zeros(xpopspace_dim + 1), xpopspace_index_range)
    result[0] = c
    result
end

"""
x^(x_power) p^(p_power)
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
            return " - $c"
        end
    end
    if real(c) == 0
        return " + $(imag(c))im"
    end
    string(c)
end

function xpopstr_stringify(coefficients::OffsetArray)
    result_terms = []
    for idx in xpopspace_index_range
        x_power = index_to_xpower(idx)
        p_power = index_to_ppower(idx)
        if coefficients[idx] != 0
            push!(result_terms, 
                "$(coefficient_format(coefficients[idx])) $(power_str("x", x_power)) $(power_str("p", p_power))")
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
end

"""
x^(x_power) * idx * p^(p_power)
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
x^(x_power) * coefficients * p^(p_power)
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
# because this function is used to calculate the normal ordering)
"""
function xpopstr_comm(x_power, p_power)
    n = x_power
    m = p_power
    result = zeros(xpopspace_dim)
    for k in 1 : min(m, n)
        numerator = - (- im)^k * factorial(n) * factorial(m)
        denominator = factorial(k) * factorial(n - k) * factorial(m - k)
        result[xpopstr_index(n - k, m - k)] = numerator / denominator
    end
    result
end

"""
Normal ordered version of x^(x_power_1) p^(p_power_1) x^(x_power_2) p^(p_power_2)
"""
function xpopstr_normal_ord(x_power_1, p_power_1, x_power_2, p_power_2)
    
end

##

@objective(model, Min, 2M[si(0), si(2)] + 3g * M[si(0), si(4)])
optimize!(model)