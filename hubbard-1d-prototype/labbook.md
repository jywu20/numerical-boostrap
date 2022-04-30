# 2022.4.30

先从一维Hubbard做起吧。

在改代码期间我强烈感到有必要写格点无关的程序，不过这些代码应该也不会用到别的地方了，so……

使用`2022-4-30-run-1.sh`提交程序。目前的配置：
```julia
const working_path = "./"
const output_name = "2022-4-30-run-1-res"

const full_output_name = working_path * output_name

# If there exists working_path * output_name already, throw an error
const no_conflict = false

# Display operators involved in the bootstrap process
const show_hubbard_opstr_basis = false 

#endregion

#region model parameter and cutoff

# The hopping parameter and local repulsion  
U = 4.0
t = 1.0

# l(O) ≤ K cutoff
K = 5
site_num = 2K + 1

# When this flag is `true`, no actual optimization will be done. For debugging only.
no_optimization = false 
max_iter = 10000
```

提示infeasible。

应注意此处$K$不能小于4，否则哈密顿量无法在算符空间中表示出来。

好像$K=4$也是会出问题的。
```
ERROR: LoadError: MethodError: Cannot `convert` an object of type Nothing to an object of type Array{Complex{Float64},1}
Closest candidates are:
  convert(::Type{Array{T,N}}, ::StaticArrays.SizedArray{S,T,N,N,Array{T,N}}) where {S, T, N} at C:\Users\wujin\.julia\packages\StaticArrays\N3rrI\src\SizedArray.jl:120
  convert(::Type{Array{T,N}}, ::StaticArrays.SizedArray{S,T,N,M,TData} where TData<:AbstractArray{T,M} where M) where {T, S, N} at C:\Users\wujin\.julia\packages\StaticArrays\N3rrI\src\SizedArray.jl:114
  convert(::Type{T}, ::AbstractArray) where T<:Array at array.jl:554
  ...
Stacktrace:
 [1] setindex!(::Array{Array{Complex{Float64},1},2}, ::Nothing, ::Int64, ::Int64) at .\array.jl:849
 [2] top-level scope at d:\Projects\numerical-boostrap\hubbard-1d-prototype\operator-algebra.jl:158        
in expression starting at d:\Projects\numerical-boostrap\hubbard-1d-prototype\operator-algebra.jl:152  
```