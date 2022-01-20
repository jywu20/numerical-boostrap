using LinearAlgebra
using Arpack
using Plots
using Flux
using Flux.Optimise: update!

##

Δx = 0.05
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

# Building the Hamiltonian 

H = - ∇² + V

##
# Variation ansatz 1

σ² = [1.0]
ψ(x) = exp(- x^2 / (2 * σ²[1]))

θ = params(σ²)

function energy()
    ψ_vec = ψ.(xs)
    ψ_vec' * H * ψ_vec / (ψ_vec' * ψ_vec)
end

grads = gradient(energy, θ)

##

opt = Descent(0.1) 
for _ in 1 : 10
    update!(opt, σ², grads[σ²])
    println(energy())
end

##
# Variation ansatz 2

σ² = [1.0]
an = [2.0, 0.0]

ψ(x) = exp(- x^2 / (2 * σ²[1])) * (an[1] + an[2] * x^2)

θ = params(σ², an)

function energy()
    ψ_vec = ψ.(xs)
    ψ_vec' * H * ψ_vec / (ψ_vec' * ψ_vec)
end

grads = gradient(energy, θ)

##

opt = Descent(0.01) 
for _ in 1 : 10
    for p in [σ², an]
        update!(opt, p, grads[p])
        println(energy())
    end
end

