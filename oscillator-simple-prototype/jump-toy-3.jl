using JuMP, CSDP
using LinearAlgebra

model = Model(CSDP.Optimizer)
@variable(model, x1)
@variable(model, x2)
@variable(model, M[1:3, 1:3], PSD)

@constraint(model, M[1, 1] == x1)
@constraint(model, M[1, 2] == x1 + im * x2)
@constraint(model, M[1, 3] == x2 + 0.5)
@constraint(model, M[2, 2] == 5x2)
@constraint(model, M[2, 3] == 0)
@constraint(model, M[3, 3] == x1 + x2)

@objective(model, Min, 2x1 + 3x2)
optimize!(model)

@show objective_value(model)
@show value(a)
@show value(b)