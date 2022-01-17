import Base.copy

#region The discrete path integral configuration of a Z2 gauge field.

"""
Z2 gauge field, **D**iscrete **P**ath **I**ntegral configuration.

`V` is the data type for spins defined on bonds. It must be something that can be multiplied like an integer, 
since 
"""
struct Z2GaugeFieldDPI{L <: AbstractLatticeWithPlaquattes, V} <: AbstractDiscretePathIntegralConfiguration
    n_τ::Int
    lattice::L
    data::Array{V, 2}
end

function copy(σ::Z2GaugeFieldDPI{L, V}) where {L <: AbstractLatticeWithPlaquattes, V}
    Z2GaugeFieldDPI{L, V}(σ.n_τ, σ.lattice, copy(σ.data))
end

field_lattice(σ::Z2GaugeFieldDPI) = σ.lattice
time_step_number(σ::Z2GaugeFieldDPI) = σ.n_τ
time_steps(σ::Z2GaugeFieldDPI) = 1 : σ.n_τ
timeslice(σ::Z2GaugeFieldDPI, τ) = σ.data[:, τ]

"""
Get value of σ^z_{ij}(τ).
"""
getindex(σ::Z2GaugeFieldDPI, i, j, τ) = σ.data[sites_to_bond(σ.lattice, i, j), τ]
"""
Set σ^z_{ij}(τ) = v.
"""
setindex!(σ::Z2GaugeFieldDPI, v, i, j, τ) = σ.data[sites_to_bond(σ.lattice, i, j), τ] = v
"""
Get value of σ^z_b(τ).
"""
getindex(σ::Z2GaugeFieldDPI, b, τ) = σ.data[b, τ]
"""
Set σ^z_b(τ) = v.
"""
setindex!(σ::Z2GaugeFieldDPI, v, b, τ) = σ.data[b, τ] = v

function ones_Z2_gauge_field_DPI(::Type{V}, lattice::L, n_τ) where {L <: AbstractLatticeWithPlaquattes, V}
    Z2GaugeFieldDPI{L, V}(n_τ, lattice, ones(V, bond_number(lattice) , n_τ))
end

function random_Z2_gauge_field_DPI(::Type{Int}, lattice::L, n_τ) where {L <: AbstractLatticeWithPlaquattes, V}
    Z2GaugeFieldDPI{L, Int}(n_τ, lattice, rand((-1, 1), bond_number(lattice), n_τ))
end

#endregion

#region The Ising gauge theory

"""
The B_p = ∏_{b ∈ ☐_p} σᶻ_b(τ) operator defined on a time step τ.

`p` is the plaquatte.
"""
prod_σ_plaquatte(σ::Z2GaugeFieldDPI, p, τ) = map(b -> σ[b, τ], plaquatte_to_bonds(field_lattice(σ), p)) |> prod

"""
∑_p B_p / N.
"""
flux_average(σ::Z2GaugeFieldDPI, τ) = mean(map(x -> prod_σ_plaquatte(σ, x, τ), sites(σ.lattice)))

function Δ_plaquatte_term(σ::Z2GaugeFieldDPI, b, τ)
    lattice = σ.lattice
    - 2 * sum(map(p -> prod_σ_plaquatte(σ, p, τ), plaquatte_containing_bond(lattice, b)))
end

function Δ_temporal_correlation_term(σ::Z2GaugeFieldDPI, b, τ)
    n_τ = time_step_number(σ)
    - 2 * σ[b, τ] * (σ[b, back_into_range(τ + 1, n_τ)] + σ[b, back_into_range(τ - 1, n_τ)])
end

"""
Parameters for discrete path integral Metropolis Monte Carlo simulation of an Ising gauge theory, i.e.
H = - J ∏_{b ∈ ☐_p} σᶻ_b(τ) - h ∏_{b ∈ ☐_p} σᶻ_b(τ).
"""
struct IsingGaugeTheoryDPIMetropolisMC{F <: AbstractFloat} <: AbstractModel
    J::F
    h::F
    Δτ::F
    β::F
    J_xy::F
    J_τ::F
end

function IsingGaugeTheoryDPIMetropolisMC(::Type{F}, σ::Z2GaugeFieldDPI, J::F, h::F, Δτ::F) where {F <: AbstractFloat}
    n_τ = time_step_number(σ)
    β = Δτ * n_τ
    J_xy = Δτ * J
    J_τ = atanh.(exp.(- 2 * Δτ * h))
    IsingGaugeTheoryDPIMetropolisMC{F}(J, h, Δτ, β, J_xy, J_τ)
end

function weight_ratio(model::IsingGaugeTheoryDPIMetropolisMC, σ::Z2GaugeFieldDPI, b::Int, τ)
    J_xy = model.J_xy
    J_τ = model.J_τ
    exp(J_xy * Δ_plaquatte_term(σ, b, τ) + J_τ * Δ_temporal_correlation_term(σ, b, τ))
end

function sweep!(model::IsingGaugeTheoryDPIMetropolisMC, σ::Z2GaugeFieldDPI, n_sweep::Integer; 
    observe = nothing, observable_type::DataType = Any)
    results = observable_type[]
    lattice_bonds = bonds(σ.lattice)

    for _ in 1 : n_sweep
        for τ in time_steps(σ)
            for b in lattice_bonds
                if rand() < weight_ratio(model, σ, b, τ)
                    σ[b, τ] *= -1 
                end
            end
        end
    end

    if observe !== nothing
        push!(results, observe(σ))
    end

    results
end

#endregion