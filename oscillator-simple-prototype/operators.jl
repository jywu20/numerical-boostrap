using SparseArrays
import Base.setindex!, Base.getindex, Base.*

"""
An operator string made of x and p. x is placed in front of p.
"""
struct XPOpString 
    # Maximal power + 1 (for example if K = 10, x^0, x^1, ..., x^9 are allowed)
    K::Int
    coefficient::SparseMatrixCSC{ComplexF64}
end

function XPOpString(K::Int)
    XPOpString(X, spzeros(K, K))
end

function XPOpString(K::Int, x_order, p_order)
    res = XPOpString(K)
    res.coefficient[x_order, p_order] = 1.0
end

function XPOpString(K::Int, x_order, p_order, val)
    res = XPOpString(K)
    res.coefficient[x_order, p_order] = val
end

function getindex(ops::XPOpString, x_order, p_order)
    ops.coefficient[x_order, p_order]
end

function setindex!(ops::XPOpString, val, x_order, p_order)
    ops.coefficient[x_order, p_order] = val
end

"""
OPE by commutation relations.
"""
function ope(a::XPOpString, b::XPOpString)
    
end