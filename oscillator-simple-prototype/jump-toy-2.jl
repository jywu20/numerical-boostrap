using JuMP, CSDP
using LinearAlgebra

model = Model(CSDP.Optimizer)
@variable(model, a)
@variable(model, b)
@variable(model, M[1:3, 1:3], PSD)

@constraint(model, M[1, 1] == a)
@constraint(model, M[1, 2] == b)
@constraint(model, M[1, 3] + a + b + 0.5 == 0)
@constraint(model, M[2, 2] == 2a + 3b + 1)
@constraint(model, M[2, 3] + 4a == 0)
@constraint(model, M[3, 3] == 4a + 5b)

@objective(model, Min, a + 2b)
optimize!(model)

@show objective_value(model)
@show value(a)
@show value(b)