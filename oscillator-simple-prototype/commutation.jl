using LinearAlgebra
using SparseArrays
import Base.setindex!, Base.getindex, Base.*, Base.+, Base.-, Base.^, Base.show, Base.==

index_to_number(x) = x - 1
number_to_index(x) = x + 1

"""
An operator string made of x and p. x is placed in front of p.
"""
struct XPOpString 
    # Maximal power + 1 (for example if K = 10, x^0, x^1, ..., x^9 are allowed)
    K::Int
    coefficient::SparseMatrixCSC{ComplexF64}
end

const default_maximal_power = 20

function XPOpString(K::Int)
    XPOpString(K, spzeros(K, K))
end

index_to_number(x) = x - 1
number_to_index(x) = x + 1

function XPOpString(coefficient::AbstractMatrix)
    K1, K2 = size(coefficient)
    if K1 != K2
        error("The coefficient matrix must be a square matrix.")
    end
    K = K1 
    XPOpString(K, sparse(coefficient))
end

function XPOpString(K::Int, x_order, p_order)
    res = XPOpString(K)
    res.coefficient[x_order, p_order] = 1.0
    res
end

function XPOpString(K::Int, x_order, p_order, val)
    res = XPOpString(K)
    res.coefficient[x_order, p_order] = val
    res
end

function getindex(ops::XPOpString, x_order, p_order)
    ops.coefficient[x_order, p_order]
end

function setindex!(ops::XPOpString, val, x_order, p_order)
    ops.coefficient[x_order, p_order] = val
end

function change_maximal_power(ops::XPOpString, K_new::Integer)
    K = ops.K
    if K == K_new
        return ops
    end
    if K > K_new
        return XPOpString(K_new, ops.coefficient[1:K_new, 1:K_new])
    end

    return XPOpString(K_new, [
        ops.coefficient          spzeros(K, K_new - K) ;
        spzeros(K_new - K, K)    spzeros(K_new - K, K_new - K)
    ])
end

==(o1::XPOpString, o2::XPOpString) = o1.coefficient == o2.coefficient

+(o1::XPOpString, o2::XPOpString) = XPOpString(o1.K, o1.coefficient + o2.coefficient)

-(o1::XPOpString, o2::XPOpString) = o1 + (-1) * o2
-(o::XPOpString) = XPOpString(o.K, - o.coefficient)

x = XPOpString(default_maximal_power, number_to_index(1), number_to_index(0))
p = XPOpString(default_maximal_power, number_to_index(0), number_to_index(1))
const_op = XPOpString(default_maximal_power, number_to_index(1), number_to_index(1))

function power_str(var, pow)
    if pow == 0
        return ""
    end
    if pow == 1
        return var 
    end
    return "$var^$pow"
end

function show(io::IO, ::MIME"text/plain", ops::XPOpString)
    terms = String[]
    row_idxs, col_idxs, entries = findnz(ops.coefficient)

    for i in 1 : length(row_idxs)
        row_idx = row_idxs[i]
        col_idx = col_idxs[i]
        entry = entries[i]
        push!(terms, "$entry $(power_str("x", index_to_number(row_idx))) $(power_str("p", index_to_number(col_idx)))")
    end
    print(io, join(terms, " + "))
end

*(x::Number, o::XPOpString) = XPOpString(o.K, x * o.coefficient)
*(o::XPOpString, x::Number) = XPOpString(o.K, x * o.coefficient)

function *(ops1::XPOpString, ops2::XPOpString)
    row_idxs1, col_idxs1, entries1 = findnz(ops1.coefficient)
    row_idxs2, col_idxs2, entries2 = findnz(ops2.coefficient)

    K = max(ops1.K, ops2.K)
    result = 0.0 * XPOpString(K, number_to_index(0), number_to_index(0))

    for i in 1 : length(row_idxs1), j in 1 : length(row_idxs2)
        x_num_1 = index_to_number(row_idxs1[i])
        p_num_1 = index_to_number(col_idxs1[i])
        x_num_2 = index_to_number(row_idxs2[j])
        p_num_2 = index_to_number(col_idxs2[j])
        coefficient1 = entries1[i]
        coefficient2 = entries2[j]

        result += XPOpString(K, number_to_index(x_num_1 + x_num_2), number_to_index(p_num_1 + p_num_2))
        for k in 1 : min(p_num_1, x_num_2)
            prefactor_num = (-im)^k * factorial(x_num_2) * factorial(p_num_1) 
            prefactor_den = factorial(k) * factorial(x_num_2 - k) * factorial(p_num_1 - k)
            prefactor = coefficient1 * coefficient2 * prefactor_num / prefactor_den
            result += prefactor * XPOpString(K, number_to_index(x_num_2 - k), number_to_index(p_num_1 - k))
        end
    end

    result
end

^(ops::XPOpString, n::Integer) = prod(repeat([ops], n))

comm(ops1::XPOpString, ops2::XPOpString) = ops1 * ops2 - ops2 * ops1

