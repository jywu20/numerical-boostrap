# Taken from https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/blob/master/test/runtests.jl

# Test cases of semidefinite programming
using NonconvexSemidefinite, NonconvexIpopt, LinearAlgebra, Test
using Distributions, DistributionsAD, ChainRulesTestUtils, Random

# Interaction strength
g = 1.0

sd_constraint((x⁴, x²)) = [
    1.0       0        x² ;
    0         x²       0  ;
    x²        0        x⁴ 
]

# Objective function
f((x⁴, x²)) = 2x² + 3 * g * x⁴ 

model = Model(f)

x0 = [0.774, 0.298]
lbs = [0.0, 0.0]
ubs = [Inf, Inf]
addvar!(model, lbs, ubs)

add_sd_constraint!(model, sd_constraint)

alg = SDPBarrierAlg(sub_alg=IpoptAlg())
options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200))
result = optimize(model, alg, x0, options = options)

E_min, minimizer, optimal_ind = result.minimum, result.minimizer, result.optimal_ind
m_mat_minimizer = sd_constraint(minimizer)
x² = minimizer[2][1]

println("ground state energy  =  $E_min")
println("⟨x²⟩                  = \n $x²")
println(eigen(m_mat_minimizer).values)
