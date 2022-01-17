# Calculate Green functions from B matrices, using a numerical stabilization scheme.

function udv_decompose(m)
    U, S, Vt = svd(m)
    (U, diagm(S), Vt')
end

function B_τ_0_udv(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, τ)
    if τ == 0
        return (I, I, I)
    end
    U, D, V = udv_decompose(B_mat(model, aux, 1))
    for τ′ in 2 : τ 
        Up, Dp, Vp = udv_decompose(B_mat(model, aux, τ′) * U * D)
        copy!(V, Vp * V)
        copy!(D, Dp)
        copy!(U, Up)
    end
    (U, D, V)
end

function B_β_τ_vdu(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, τ)
    if τ == n_τ
        return (I, I, I)
    end
    V, D, U = udv_decompose(B_mat(model, aux, n_τ))
    for τ′ in n_τ - 1 : -1 : τ + 1
        Vp, Dp, Up = udv_decompose(D * U * B_mat(model, aux, τ′))
        copy!(V, V * Vp)
        copy!(D, Dp)
        copy!(U, Up)
    end
    (V, D, U)
end

function G_τ_τ(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, τ)
    Ur, Dr, Vr = B_τ_0_udv(model, aux, τ)
    Vl, Dl, Ul = B_β_τ_vdu(model, aux, τ)
    U, D, V = udv_decompose(inv(Ul * Ur) + Dr * (Vr * Vl) * Dl)
    inv(V * Ul) * inv(D) * inv(Ur * U)
end

function G_τ_τ!(model::Z2SpinlessFermionSimpleDQMC, aux::Z2SpinlessFermionSimpleAuxField, τ)
    aux.G[:, :, τ] = G_τ_τ(model, aux, τ)
end