using LinearAlgebra
using Convex, SCS

##

# The total energy
E = Variable()
# ⟨x²⟩
x² = Variable()

# Interaction strength
g = 1.0

# Maximal K in (7) in 2004.10212
K = 7

# Recursively calculate ⟨xⁿ⟩

function expected_xn(n)
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
        return (E - 2x²) / 3g # According to (3) in 2004.10212    
    end
    
    t = n - 3
    # According to (6) in 2004.10212
    (4t * E * expected_xn(t - 1) + t * (t - 1) * (t - 2) * expected_xn(t - 3) - 4 * (t + 1) * expected_xn(t + 1)) / (4g * (t + 2))
end

# Construct the optimalization problem 

# Here we are unable to use the expression
# M = [expected_xn(i + j) for i in 0 : K, j in 0 : K]
# directly because Convex.jl does not recognize a matrix of variables as a variable.
# We need a conversion process.

function convert_to_variable(arr::AbstractMatrix{<:Convex.AbstractExprOrValue})
    x = Variable(size(arr))
    for i = 1:size(arr, 1), j = 1:size(arr, 2)
        add_constraint!(x, x[i,j] == arr[i,j])
    end
    return x
end

M = convert(Matrix{Convex.AbstractExprOrValue}, [expected_xn(i + j) for i in 0 : K, j in 0 : K])
M = convert_to_variable(M)

constraint = (M in :SDP)
problem = minimize(E, constraint)

##
# Solving the problem.

solve!(problem, SCS.Optimizer())

##
