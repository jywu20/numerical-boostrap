using JuMP, CSDP
using LinearAlgebra

model = Model(CSDP.Optimizer)
@variable(model, x1)
@variable(model, x2)
@variable(model, M[1:6, 1:6], PSD)

I22 = Matrix(1*I, 2, 2)
Im22 = [
    0   -1 ; 
    1    0
]

@constraint(model, M[1:2, 1:2] .== x1 * I22)
@constraint(model, M[1:2, 3:4] .== x1 * I22 + x2 * Im22)
@constraint(model, M[1:2, 5:6] .== (x2 + 0.5) * I22)
@constraint(model, M[3:4, 3:4] .== 5x2 * I22)
@constraint(model, M[3:4, 5:6] .== 0 * I22)
@constraint(model, M[5:6, 5:6] .== (x1 + x2) * I22)

@objective(model, Min, 2x1 + 3x2)
optimize!(model)

@show objective_value(model)
@show value(x1)
@show value(x2)