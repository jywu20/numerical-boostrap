using JuMP, CSDP
using LinearAlgebra

model = Model(CSDP.Optimizer)
@variable(model, a)
@variable(model, b)
@variable(model, c)

M = [a b;
     b c]

@SDconstraint(model, M ⪰ zeros(2, 2))
@constraint(model, a + b == c)

@objective(model, Min, a + 2b + c)
optimize!(model)