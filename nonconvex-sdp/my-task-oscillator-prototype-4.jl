# Taken from https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/blob/master/test/runtests.jl

# Test cases of semidefinite programming
using NonconvexSemidefinite, NonconvexIpopt, LinearAlgebra, Test
using Distributions, DistributionsAD, ChainRulesTestUtils, Random

# Interaction strength
g = 1.0
# Maximal K in (7) in 2004.10212
K = 7

# Objective function
f(x) = 2x[2] + 3 * g * x[4]
model = Model(f)

function expected_xn_def(x⁴, x², n)
    if isodd(n)
        return 0.0
    end
    if n == 0
        return 1.0
    end
    if n == 2
        return x²
    end
    if n == 4
        return x⁴  
    end
    
    t = n - 3
    # According to (6) in 2004.10212
    (4t * (2x² + 3 * g * x⁴) * expected_xn_def(x⁴, x², t - 1) + t * (t - 1) * (t - 2) * expected_xn_def(x⁴, x², t - 3) - 4 * (t + 1) * expected_xn_def(x⁴, x², t + 1)) / (4g * (t + 2))
end

function power_to_index(n)
    Int(n / 2)
end

function expected_xn(x, n)
    if isodd(n)
        return 0.0
    end
    if n == 0
        return 1.0
    end
    x[power_to_index(n)]
end

# The initial value of E and ⟨x²⟩ are [1.37, 0.298], and therefore 
# x⁴ = (E - 2x²) / 3g = 0.258 
x²_0 = 0.298
x⁴_0 = 0.258
x0 = zeros(K)
x0[power_to_index(2)] = x²_0
x0[power_to_index(4)] = x⁴_0
for n in 6 : 2 : 2K
    x0[power_to_index(n)] = expected_xn_def(x⁴_0, x²_0, n)
end

# Register K variables
lbs = zeros(K)
ubs = fill(Inf, K)
addvar!(model, lbs, ubs)

n = K
while true
    t = n - 3
    if t < 0
        break
    end
    # According to (6) in 2004.10212
    constraint_func = (x) -> (
        4t * (2x[power_to_index(2)] + 3 * g * x[power_to_index(4)]) * expected_xn(x, t - 1) 
        + t * (t - 1) * (t - 2) * expected_xn(x, t - 3) 
        - 4 * (t + 1) * expected_xn(x, t + 1)
        - 4g * (t + 2) * expected_xn(x, t + 3)
    )
    add_eq_constraint!(model, constraint_func)
    n -= 2
end

sd_constraint(x) = [expected_xn(x, i + j) for i in 0 : K, j in 0 : K]
add_sd_constraint!(model, sd_constraint)

alg = SDPBarrierAlg(sub_alg=IpoptAlg())
options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200))
result = optimize(model, alg, x0, options = options)

E_min, minimizer, optimal_ind = result.minimum, result.minimizer, result.optimal_ind
m_mat_minimizer = sd_constraint(minimizer)
x² = minimizer[power_to_index(2)]

println("ground state energy  =  $E_min")
println("⟨x²⟩                  = \n $x²")
