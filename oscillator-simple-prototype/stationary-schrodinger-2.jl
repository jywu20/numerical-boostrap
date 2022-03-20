using LinearAlgebra
using Arpack
using Plots
using OffsetArrays
using JLD2

Δx = 0.02
L = 10
# The list of all x coordinate points included 
xs = - L : Δx : L
xs_len = length(xs)

# The Laplacian operator. The boundary condition is integrated into the matrix representation.
# p² = - ∇²
∇² = zeros(xs_len, xs_len)
∇²[1, 1] = -2
∇²[1, 2] = 1
for i in 2 : xs_len - 1
    ∇²[i, i - 1] = 1
    ∇²[i, i] = -2
    ∇²[i, i + 1] = 1
end
∇²[xs_len, xs_len] = -2
∇²[xs_len, xs_len - 1] = 1
∇² /= Δx^2

# The potential
g = 1
V = diagm(xs.^2 + g * xs.^4)

# Building the Hamiltonian and diagonalize

H = - ∇² + V
result = eigs(H, nev = 10, which = :SM)

##

ψ_ground = result[2][:, 1]
plot(xs, ψ_ground, legend = false)

##

L_max = 5

∇ = zeros(xs_len, xs_len)
for i in 1 : xs_len - 1
    ∇[i, i] = -1
    ∇[i, i + 1] = 1
end
∇[xs_len, xs_len] = -1
∇ /= Δx
p_op = - im * ∇

# First index: x_power, second index: p_power
xpopstr_expected_ode = OffsetArray(zeros(ComplexF64, 2L_max + 1, 2L_max + 1), 0 : 2L_max, 0 : 2L_max)
x_op = diagm(xs) 

for x_power in 0 : 2L_max
    for p_power in 0 : 2L_max
        xpopstr_expected_ode[x_power, p_power] = 
            ψ_ground' * x_op^x_power * p_op^p_power * ψ_ground / (ψ_ground' * ψ_ground)
    end
end

##

working_path = "D:\\Projects\\numerical-boostrap\\oscillator-simple-prototype\\"
@save working_path * "xpopstr_expected_ode-dx-0.02-l-10.jld2" xpopstr_expected_ode