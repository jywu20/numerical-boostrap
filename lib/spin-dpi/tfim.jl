#region The discrete path integral configuration of an Ising field.

struct IsingFieldDPI{L <: AbstractLattice, V} <: AbstractDiscretePathIntegralConfiguration
    n_τ::Int
    lattice::L
    data::Array{V, 2}
end

field_lattice(s::IsingFieldDPI) = s.lattice
time_step_number(s::IsingFieldDPI) = s.n_τ
time_steps(s::IsingFieldDPI) = 1 : s.n_τ
timeslice(s::IsingFieldDPI, τ) = s.data[:, τ]

setindex!(s::IsingFieldDPI, v, i, τ) = s.data[i, τ] = v
getindex(s::IsingFieldDPI, i, τ) = s.data[i, τ]

function ones_Ising_field_DPI(::Type{V}, lattice::L, n_τ) where {L <: AbstractLattice{Int}, V}
    IsingFieldDPI{L, V}(n_τ, lattice, ones(V, site_number(lattice) , n_τ))
end

#endregion

#region Discrete path integral Monte Carlo for transverse field Ising model

function magnetization(s::IsingFieldDPI, τ)
    abs(sum(s.data[:, τ])) / site_number(s.lattice)
end

function Δ_spacial_correlation_term(s::IsingFieldDPI, i, τ)
    lattice = field_lattice(s)
    nearest_site_values = map(j -> s[j, τ], nearest_neighbors(lattice, i))
    - 2 * s[i, τ] * sum(nearest_site_values)
end

function Δ_temporal_correlation_term(s::IsingFieldDPI, i, τ)
    n_τ = time_step_number(s)
    - 2 * s[i, τ] * (s[i, back_into_range(τ + 1, n_τ)] + s[i, back_into_range(τ - 1, n_τ)])
end

"""
Parameters for discrete path integral Metropolis Monte Carlo simulation of the transverse field Ising model
H = - J ∑_{⟨i,j⟩} sᶻ_i sᶻ_j - h ∑_i sˣ_i.
"""
struct TransverseFieldIsingModelDPIMetropolisMC{F <: AbstractFloat} <: AbstractModel
    J::F
    h::F
    Δτ::F
    β::F
    J_xy::F
    J_τ::F
end

function TransverseFieldIsingModelDPIMetropolisMC(::Type{F}, s::IsingFieldDPI, J::F, h::F, Δτ::F) where {F <: AbstractFloat}
    n_τ = time_step_number(s)
    β = Δτ * n_τ
    J_xy = Δτ * J
    J_τ = atanh.(exp.(- 2 * Δτ * h))
    TransverseFieldIsingModelDPIMetropolisMC{F}(J, h, Δτ, β, J_xy, J_τ)
end

function weight_ratio(model::TransverseFieldIsingModelDPIMetropolisMC, s::IsingFieldDPI, i, τ)
    J_xy = model.J_xy
    J_τ = model.J_τ
    exp(J_xy * Δ_spacial_correlation_term(s, i, τ) + J_τ * Δ_temporal_correlation_term(s, i, τ))
end

"""
The first parameter is only relevant for its type.
"""
function update!(_::TransverseFieldIsingModelDPIMetropolisMC, s::IsingFieldDPI, i, τ)
    s[i, τ] *= -1
end

function sweep!(model::TransverseFieldIsingModelDPIMetropolisMC, s::IsingFieldDPI, n_sweep::Integer; 
    observe = nothing, observable_type::DataType = Any)
    results = observable_type[]
    lattice_sites = sites(s.lattice)

    for _ in 1 : n_sweep
        for τ in time_steps(s)
            for i in lattice_sites
                if rand() < weight_ratio(model, s, i, τ)
                    s[i, τ] *= -1 
                end
            end
        end
    end

    if observe !== nothing
        push!(results, observe(s))
    end

    results
end

#endregion