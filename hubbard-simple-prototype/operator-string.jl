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
# Whether there is no operator, a spin-up operator, a spin-down operator or both a spin-up and a spin-down 
# operator on each site. The index of current_op_str is the same as the index of site_list
current_op_str = [:up]
# For an operator string with r operators, its size is r + ∑_i |x_i|_1, where x_i is the coordinate of 
# the i-th operator
current_op_str_size = 1
current_site = 1

# To avoid out of memory error, we set a max qualified operator string limit.
max_qualified_opstr_num = 10000
# Tp avoid running for too long
max_search_num = 10000

qualified_opstr = Vector{Symbol}[]

finished = false
search_count = 0
while ! finished 

    global search_count += 1

    if length(qualified_opstr) > max_qualified_opstr_num || search_count > max_search_num
        @warn "max_qualified_opstr_num reached."
        break
    end
    
    if current_site < site_num && current_op_str_size ≤ K
        println("$current_site, add site, $current_op_str")
        # Add a site
        push!(current_op_str, :no)
        current_site += 1
        # No increasing in current_op_str_size, because no electron operatore is actually added
        continue
    end

    if current_site == site_num && current_op_str_size ≤ K
        println("$current_site, record, $current_op_str")
        push!(qualified_opstr, copy(current_op_str))
    end

    # Move to the next configuration: from no operator to spin-up operator
    if current_op_str[current_site] == :no
        println("$current_site, no to up, $current_op_str")
        current_op_str[current_site] = :up 
        # Operator num + 1, and the 1-norm of the current site must also be added 
        current_op_str_size += 1 + site_norm_1_from_center[current_site]           
        continue
    end

    # Move to the next configuration: from no operator to spin-up operator
    if current_op_str[current_site] == :up
        println("$current_site, up to dn, $current_op_str")
        # No change in the coordinates and number of operators, so no change on current_op_str_size
        current_op_str[current_site] = :dn
        continue
    end

    # Move to the next configuration: from spin-down operator to both spin up and spin down
    if current_op_str[current_site] == :dn
        println("$current_site, dn to both, $current_op_str")
        current_op_str[current_site] = :both
        # Operator num + 1, and the 1-norm of the current site must also be added 
        current_op_str_size += 1 + site_norm_1_from_center[current_site]
        continue
    end

    if current_op_str[current_site] == :both

        println("$current_site, moving back, $current_op_str")
        
        while current_op_str[current_site] == :both
            current_op_str_size -= 2 + 2site_norm_1_from_center[current_site]
            current_site -= 1
            pop!(current_op_str)

            if current_site == 0
                global finished = true
                break
            end
        end

        if finished
            break
        end

        # Move to the next configuration: from no operator to spin-up operator
        if current_op_str[current_site] == :no
            println("$current_site, after moving back no to up, $current_op_str")
            current_op_str[current_site] = :up 
            # Operator num + 1, and the 1-norm of the current site must also be added 
            current_op_str_size += 1 + site_norm_1_from_center[current_site]           
            continue
        end

        # Move to the next configuration: from no operator to spin-up operator
        if current_op_str[current_site] == :up
            println("$current_site, after moving back up to dn, $current_op_str")
            # No change in the coordinates and number of operators, so no change on current_op_str_size
            current_op_str[current_site] = :dn
            continue
        end

        # Move to the next configuration: from spin-down operator to both spin up and spin down
        if current_op_str[current_site] == :dn
            println("$current_site, after moving back dn to both, $current_op_str")
            current_op_str[current_site] = :both
            # Operator num + 1, and the 1-norm of the current site must also be added 
            current_op_str_size += 1 + site_norm_1_from_center[current_site]
            continue
        end
    end
end