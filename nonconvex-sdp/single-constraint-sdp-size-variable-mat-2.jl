# Taken from https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/blob/master/test/runtests.jl

# Test cases of semidefinite programming
using NonconvexSemidefinite, NonconvexIpopt, LinearAlgebra, Test
using Distributions, DistributionsAD, ChainRulesTestUtils, Random

K = 3
g = 1

Mij(x⁴, x², n) = begin
    if n == 0
        return 1.0
    end
    if isodd(n)
        return 0.0
    end
    if n == 2
        return x²
    end
    if n == 4
        return x⁴
    end
    
    (n - 4) * Mij(x⁴, x², n - 2) * x² + Mij(x⁴, x², n - 4)
end

sd_constraint((x⁴, x²)) = [Mij(x⁴, x², i + j) for i in 0 : K, j in 0 : K]

# Objective function
f((x⁴, x²)) = 2x² + 3x⁴

model = Model(f)

x0 = [1.0, 0.8]
lbs = [0.0, 0.0]
ubs = [Inf, Inf]
addvar!(model, lbs, ubs)

add_sd_constraint!(model, sd_constraint)

alg = SDPBarrierAlg(sub_alg=IpoptAlg())
options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=400))
result = optimize(model, alg, x0, options = options)

minimum, minimizer, optimal_ind = result.minimum, result.minimizer, result.optimal_ind
m_mat_minimizer = sd_constraint(minimizer)

@show minimum
@show eigen(m_mat_minimizer).values
