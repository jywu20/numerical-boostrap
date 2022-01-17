# Z2 - Ising coupling

struct Z2IsingMinimalCouplingDPI{F <: AbstractFloat}
    J::F
    Δτ::F
end

"""
Flipping the Z2 field
"""
function weight_ratio_σ(model::Z2IsingMinimalCouplingDPI, σ::Z2GaugeFieldDPI, s::IsingFieldDPI, b, τ)
    lattice = σ.lattice
    J = model.J
    i, j = bond_to_sites(lattice, b)
    exp(2 * Δτ * J * (s[i, τ] + s[j, τ]) * σ[i, j, τ])
end

"""
Flipping the Ising field
"""
function weight_ratio_s(model::Z2IsingMinimalCouplingDPI, σ::Z2GaugeFieldDPI, s::IsingFieldDPI, i, τ)
    lattice = σ.lattice
    J = model.J
    Δτ = model.Δτ
    exp(2 * Δτ * J * sum(map(j -> s[j, τ] * σ[i, j, τ] , nearest_neighbors(lattice, i))))
end