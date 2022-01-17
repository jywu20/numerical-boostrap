function l1norm(coord::Tuple{Integer, Integer})
    abs(coord[1]) + abs(coord[2])
end

"""
(x, y) → (x, -y)
"""
function reflect_x(coord::Tuple{Integer, Integer})
    (coord[1], - coord[2])
end

"""
(x, y) → (-x, y)
"""
function reflect_y(coord::Tuple{Integer, Integer})
    (- coord[1], coord[2])
end

"""
(x, y) → (-y, x)
"""
function anticlockwise_90(coord::Tuple{Integer, Integer})
    (- coord[2], coord[1])
end