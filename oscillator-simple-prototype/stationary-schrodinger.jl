using LinearAlgebra
using Arpack
using Plots

##
# Space discretization 

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

# Building the Hamiltonian and diagonalize

H = - ∇² + V
result = eigs(H, nev = 10, which = :SM)

##

scatter(result[1], legend = false)

##

plot(xs, result[2][:, 1], legend = false)

## 

function spectrum_anham(L, Δx, g)
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
    V = diagm(xs.^2 + g * xs.^4)

    # Building the Hamiltonian and diagonalize

    H = - ∇² + V
    H, eigs(H, nev = 10, which = :SM), xs
end

##

L = 10
g = 1

Δx_range = [0.1, 0.05, 0.01]
E0 = []

for Δx in Δx_range
    _, result, _ = spectrum_anham(L, Δx, g)
    push!(E0, result[1][1])
end

plot(E0)

##

L = 10
g = 1

_, result, xs = spectrum_anham(L, Δx, g)
E0 = result[1][1]
ψ = result[2][:, 1]

xⁿ_expected = []

for n in 2 : 2 : 20
    push!(xⁿ_expected, ψ' * (xs.^n .* ψ))
end

xⁿ_expected