# Taken from testset "Single semi-definite constraint"
# from https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/blob/master/test/runtests.jl

# Test cases of semidefinite programming
using NonconvexSemidefinite, NonconvexIpopt, LinearAlgebra, Test
using Distributions, DistributionsAD, ChainRulesTestUtils, Random
Random.seed!(1)

# Test setting
mat_dim = 2
mat_length = mat_dim^2
n_sample = 300

function random_psd_mat(mat_dim)
    _mat = randn(mat_dim, mat_dim)
    return _mat' * _mat
end

function sd_constraint((x_L, x_D))
    decompress_symmetric(x_L, x_D)
end

# Randomize groundtruth
μ, Σ = randn(mat_dim), random_psd_mat(mat_dim)

# Generate samples
samples = rand(MvNormal(μ, Σ), n_sample)

# Objective function
function f((x_L, x_D))
    -loglikelihood(MvNormal(μ, decompress_symmetric(x_L, x_D)), samples)
end

model = Model(f)

mat_x0 = random_psd_mat(mat_dim)
x0 = [mat_x0[NonconvexSemidefinite.lowertriangind(mat_x0)], diag(mat_x0)]
lbs = [fill(-Inf, length(x0[1])), zeros(length(x0[2]))]
ubs = [fill(Inf, length(x0[1])), fill(Inf, length(x0[2]))]
addvar!(model, lbs, ubs)

add_sd_constraint!(model, sd_constraint)

alg = SDPBarrierAlg(sub_alg=IpoptAlg())
options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200))
result = optimize(model, alg, x0, options = options)

minimum, minimizer, optimal_ind = result.minimum, result.minimizer, result.optimal_ind
_Σ = sd_constraint(minimizer)

println("result: \n $result")

println("minimum: \n $minimum")
println("minimizer: \n $minimizer")
println("_Σ: \n $(_Σ)")

println("Σ: \n $Σ")
println("abs(_Σ - Σ): \n $(abs.(_Σ - Σ))")
println("mean(abs(_Σ - Σ)): \n $(mean(abs.(_Σ - Σ)))")

@test Σ ≈ _Σ rtol = 0.1
