# Taken from https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/blob/master/test/runtests.jl

# Test cases of semidefinite programming
using NonconvexSemidefinite, NonconvexIpopt, LinearAlgebra, Test
using Distributions, DistributionsAD, ChainRulesTestUtils, Random

# Interaction strength
g = 1.0
# Maximal K in (7) in 2004.10212
K = 7

function expected_xn(vars, n)
    E = vars[1][1]
    x² = vars[2][1]
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
        return (E - 2x²) / 3g # According to (3) in 2004.10212    
    end
    
    t = n - 3
    # According to (6) in 2004.10212
    (4t * E * expected_xn(vars, t - 1) + t * (t - 1) * (t - 2) * expected_xn(vars, t - 3) - 4 * (t + 1) * expected_xn(vars, t + 1)) / (4g * (t + 2))
end

sd_constraint(vars) = [expected_xn(vars, i + j) for i in 0 : K, j in 0 : K]

# Objective function
f(vars) = vars[1][1]

model = Model(f)

x0 = [[1.37], [0.298]]
lbs = [[0.0], [0.0]]
ubs = [[Inf], [Inf]]
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
