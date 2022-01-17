function Δ_mat(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, b, τ)
    σ = aux.σ
    Δ_mat_def(model, σ, b, τ)
end

function weight_ratio(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, b, τ)
    det(I + Δ_mat(model, aux, b, τ) * (I - aux.G[:, :, τ]))
end

"""
Note that after the update, equal time Green functions at other time steps are outdated.
"""
function update!(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, b, τ)
    σ = aux.σ

    Gτ = aux.G[:, :, τ]
    Δ = Δ_mat(model, aux, b, τ)
    Gτ′ = Gτ * inv(I + Δ * (I - Gτ))
    aux.G[:, :, τ] = Gτ′

    σ[b, τ] *= -1
end

"""
Propagate from τ to τ + 1
"""
function propagate_forward!(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, τ)
    if aux.propagate_count == model.n_wrap
        G_τ_τ!(model, aux, τ + 1)
        aux.propagate_count = 0
        return
    end
    Bτp1 = B_mat(model, aux, τ + 1)
    aux.G[:, :, τ + 1] = Bτp1 * aux.G[:, :, τ] * inv(Bτp1)
    aux.propagate_count += 1
end

"""
Propagate from τ to τ - 1
"""
function propagate_backward!(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, τ)
    if aux.propagate_count == model.n_wrap
        G_τ_τ!(model, aux, τ - 1)
        aux.propagate_count = 0
        return
    end
    Bτ = B_mat(model, aux, τ)
    aux.G[:, :, τ - 1] = inv(Bτ) * aux.G[:, :, τ] * Bτ
    aux.propagate_count += 1
end
