using JuMP, Ipopt, LinearAlgebra

# Interaction strength
g = 1.0
# Maximal K in (7) in 2004.10212
K = 7

function expected_xn_def(x⁴, x², n)
    if isodd(n)
        return 0.0
    end
    if n == 0
        return 1.0
    end
    if n == 2
        return x²
    end
    if n == 4
        return x⁴  
    end
    
    t = n - 3
    # According to (6) in 2004.10212
    (4t * (2x² + 3 * g * x⁴) * expected_xn_def(x⁴, x², t - 1) + t * (t - 1) * (t - 2) * expected_xn_def(x⁴, x², t - 3) - 4 * (t + 1) * expected_xn_def(x⁴, x², t + 1)) / (4g * (t + 2))
end

function power_to_index(n)
    Int(n / 2)
end

function expected_xn(x, n)
    if isodd(n)
        return 0.0
    end
    if n == 0
        return 1.0
    end
    x[power_to_index(n)]
end

problem_index_to_storage_index(i) = i + 1
si = problem_index_to_storage_index
storage_index_to_problem_index(i) = i - 1
pi = storage_index_to_problem_index

model = Model(Ipopt.Optimizer)
@variable(model, M[1 : K + 1, 1 : K + 1], PSD)

# Building the M matrix

# When i + j is odd, the corresponding matrix element is zero
for i in 0 : K
    for j in 0 : K
        if isodd(i + j) && i <= j
            @constraint(model, M[si(i), si(j)] == 0) 
        end
    end
end

@constraint(model, M[si(0), si(0)] == 0)

for n in 1 : K  
    iplusj = 2n   # Iterate over possible i+j values
    halfpoint = Int(floor(n / 2))
    for i in 1 : halfpoint
        j = n - i
        @constraint(model, M[si(0), si(n)] == M[si(i), si(j)]) 
    end
end

# If we take K = 7,
# The final result has 36 variables and 29 constraints, so there are 36 - 29 = 7 free variables, 
# the number of which is K

# Apply the poor man's OPE

n = K
while true
    t = n - 3
    if t - 3 < 0
        break
    end
    # According to (6) in 2004.10212
    @NLconstraint(model, 4t * (2M[si(0), si(2)] + 3 * g * M[si(0), si(4)]) * M[si(0), si(t - 1)] 
    + t * (t - 1) * (t - 2) * M[si(0), si(t - 3)] 
    - 4 * (t + 1) * M[si(0), si(t + 1)]
    - 4g * (t + 2) * M[si(0), si(t + 3)] == 0)
    n -= 2
end

##
# Boostrap

@objective(model, Min, 2M[si(0), si(2)] + 3g * M[si(0), si(4)])
optimize!(model)