using JuMP, CSDP, LinearAlgebra, Test

# Interaction strength
g = 1.0
# Maximal m and n in x^m x^n. The maximal length of an operator is 2K
K = 7

xpopspace_dim = (K + 1)^2

##
# Define variables. We label the expectation values as the follows: x^m p^n is labeled as m * (K + 1) + n + 1
xpopstr_index(x_power, p_power) = x_power * (K + 1) + p_power + 1
index_to_xpower(idx) = Int(floor((idx - 1) / (K + 1)))
index_to_ppower(idx) = (idx - 1) % (K + 1)

@test [xpopstr_index(i, j) for i in 0 : K for j in 0 : K] == 1 : 64

##

# Calculate [x^n, p^m] according to the McCoy's formula
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

model = Model(CSDP.Optimizer)
@variable(model, ⟨xpopstr⟩[1 : ], PSD)

# Building the M matrix

##

@objective(model, Min, 2M[si(0), si(2)] + 3g * M[si(0), si(4)])
optimize!(model)