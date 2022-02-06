# Taken from https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/blob/master/test/runtests.jl

# Test cases of semidefinite programming
using NonconvexSemidefinite, NonconvexIpopt, LinearAlgebra, Test
using Distributions, DistributionsAD, ChainRulesTestUtils, Random

# Interaction strength
g = 1.0
# Maximal K in (7) in 2004.10212
K = 7

function expected_xn(x⁴, x², n)
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
    (4t * (2x² + 3 * g * x⁴) * expected_xn(x⁴, x², t - 1) + t * (t - 1) * (t - 2) * expected_xn(x⁴, x², t - 3) - 4 * (t + 1) * expected_xn(x⁴, x², t + 1)) / (4g * (t + 2))
end

sd_constraint((x⁴, x²)) = [expected_xn(x⁴, x², i + j) for i in 0 : K, j in 0 : K]

# Objective function
f((x⁴, x²)) = 2x² + 3 * g * x⁴

model = Model(f)

# 原本是
# x0 = [1.37, 0.298]
# 在改用x²和x⁴后，初始值为
# x⁴ = (E - 2x²) / 3g = 0.258 
x0 = [0.258, 0.298]
lbs = [0.0, 0.0]
ubs = [Inf, Inf]
addvar!(model, lbs, ubs)

add_sd_constraint!(model, sd_constraint)

alg = SDPBarrierAlg(sub_alg=IpoptAlg())
options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200))
result = optimize(model, alg, x0, options = options)

E_min, minimizer, optimal_ind = result.minimum, result.minimizer, result.optimal_ind
m_mat_minimizer = sd_constraint(minimizer)
x² = minimizer[2]

println("ground state energy  =  $E_min")
println("⟨x²⟩                  = \n $x²")
