using JuMP, CSDP, LinearAlgebra, Test

# Interaction strength
g = 1.0
# Maximal m and n in x^m x^n. The maximal length of an operator is 2K
K = 7

# The -1 term comes from the fact that ⟨x^0 p^0⟩ = 1 is always the case and is not an optimization variable.
xpopspace_dim = (2K + 1)^2 - 1

##
# Define variables. We label the expectation values as the follows: x^m p^n is labeled as m * (K + 1) + n + 1
xpopstr_index(x_power, p_power) = x_power * (2K + 1) + p_power
index_to_xpower(idx) = Int(floor(idx / (2K + 1)))
index_to_ppower(idx) = idx % (2K + 1)

@test [xpopstr_index(i, j) for i in 0 : 2K for j in 0 : 2K] == 0 : xpopspace_dim

##

model = Model(CSDP.Optimizer)
@variable(model, xpopstr_expected[1 : xpopspace_dim])

# Calculate coefficients of [x^n, p^m] according to the McCoy's formula (not the operator, 
# because this function is used to calculate the normal ordering)
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

##

@objective(model, Min, 2M[si(0), si(2)] + 3g * M[si(0), si(4)])
optimize!(model)