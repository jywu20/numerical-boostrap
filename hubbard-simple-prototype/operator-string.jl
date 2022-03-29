using LinearAlgebra
using ProgressMeter

# l(O) ≤ K
K = 3 
site_num = (2K + 1)^2

site_list = Vector{Int}[]
current_site = [0, 0]
for n in 1 : K
    for i in 1 : 2n - 1
        push!(site_list, current_site)
        current_site += [0, 1]
    end
    for i in 1 : 2n - 1
        push!(site_list, current_site)
        current_site -= [1, 0]
    end
    for i in 1 : 2n
        push!(site_list, current_site)
        current_site -= [0, 1]
    end
    for i in 1 : 2n
        push!(site_list, current_site)
        current_site += [1, 0]
    end
end

for i in 1 : 2K + 1
    push!(site_list, current_site)
    current_site += [0, 1]
end

site_norm_1_from_center = map(x -> norm(x, 1), site_list)

# The set of labels of electron opeartors is site_list ⊗ spin. Whether an operator is a creation or 
# annihilation operator is not represented in the label, because we are only going to deal with 
# *normal ordered* operator strings, and therefore what are creation and annihilation operators 
# are defined for operator strings. 
spin_status_on_one_site = [:no, :up, :dn, :both]

# DFS of possible operator string
current_op_str = [:up]
current_op_str_size = 1
current_site = 1

qualified_opstr = Vector{Symbol}[]

finished = false
while ! finished 
    if current_site < site_num && current_op_str_size ≤ K
        # Add a site
        push!(current_op_str, :no)
        current_site += 1
        # No increasing in current_op_str_size, because no electron operatore is actually added
        continue
    end

    # current_site == site_num
    if current_op_str_size ≤ K
        push!(qualified_opstr, current_op_str)
    end

    # Move to the next configuration
    if current_op_str[current_site] == :no
        current_op_str[current_site] = :up 
        current_op_str_size += site_norm_1_from_center[current_site]           
        continue
    end

    if current_op_str[current_site] == :up
        current_op_str[current_site] = :dn
        continue
    end

    if current_op_str[current_site] == :dn
        current_op_str[current_site] = :both
        current_op_str_size += site_norm_1_from_center[current_site]
        continue
    end

    if current_op_str[current_site] == :both
        while current_op_str[current_site] == :both
            current_site -= 1
            current_op_str_size -= 2site_norm_1_from_center[current_site]

            if current_site == 0
                break
            end
        end
        continue
    end
end