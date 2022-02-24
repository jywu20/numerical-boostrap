Reading note of Xizhi Han's code
======

The paper mentioned is 2004.10212

# `demo.py`: the entrance

Line 22-24: building the Hamiltonian, operator algebra, etc.

Line 35: [`minimize`](#minimize) is defined from line 159 in [`optimize.py`](#optimizepy).

# `optimize.py`

The algorithm is a trust-region sequential semidefinite programming with regularization on l2 norm of the parameters.

## `minimize`

Starting from line 159.

Parameters:
- `quad_cons`: see the discussion around (14).
- `op_cons`:  