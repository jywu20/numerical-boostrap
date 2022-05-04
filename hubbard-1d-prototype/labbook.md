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
看起来这里头的问题是，用于span出M矩阵的算符中选两个做乘积得到的某些算符并不在算符空间中。

# 2022.5.4

将$K=5$时的约束和导致这些约束的算符打印出来：（见`2022-5-3-run-1-res`）
```
1     0 = 0
c†(11)     -c†(21) - c†(31) - 4 c†(1-1) c†(11) c(1-1) = 0
c†(1-1)     -c†(2-1) - c†(3-1) + 4 c†(1-1) c†(11) c(11) = 0
c†(1-1) c†(11)     4 c†(1-1) c†(11) - c†(1-1) c†(21) - c†(1-1) c†(31) + c†(11) c†(2-1) + c†(11) c†(3-1) = 0
c†(31) c(31)     -c†(11) c(31) + c†(31) c(11) + c†(31) c(51) - c†(51) c(31) = 0
c†(3-1) c(31)     -c†(1-1) c(31) + c†(3-1) c(11) + c†(3-1) c(51) - c†(5-1) c(31) = 0
c†(31) c(3-1)     -c†(11) c(3-1) + c†(31) c(1-1) + c†(31) c(5-1) - c†(51) c(3-1) = 0
c†(3-1) c(3-1)     -c†(1-1) c(3-1) + c†(3-1) c(1-1) + c†(3-1) c(5-1) - c†(5-1) c(3-1) = 0
c†(21) c(21)     -c†(11) c(21) + c†(21) c(11) + c†(21) c(41) - c†(41) c(21) = 0
c†(2-1) c(21)     -c†(1-1) c(21) + c†(2-1) c(11) + c†(2-1) c(41) - c†(4-1) c(21) = 0
c†(21) c(2-1)     -c†(11) c(2-1) + c†(21) c(1-1) + c†(21) c(4-1) - c†(41) c(2-1) = 0
c†(2-1) c(2-1)     -c†(1-1) c(2-1) + c†(2-1) c(1-1) + c†(2-1) c(4-1) - c†(4-1) c(2-1) = 0
c(11)     c(21) + c(31) - 4 c†(1-1) c(1-1) c(11) = 0
c†(11) c(11)     c†(11) c(21) + c†(11) c(31) - c†(21) c(11) - c†(31) c(11) = 0
c†(1-1) c(11)     c†(1-1) c(21) + c†(1-1) c(31) - c†(2-1) c(11) - c†(3-1) c(11) = 0
c†(1-1) c†(11) c(11)     4 c†(1-1) c†(11) c(11) + c†(1-1) c†(11) c(21) + c†(1-1) c†(11) c(31) - c†(1-1) c†(21) c(11) - c†(1-1) c†(31) c(11) + c†(11) c†(2-1) c(11) + c†(11) c†(3-1) c(11) = 0
c(1-1)     c(2-1) + c(3-1) + 4 c†(11) c(1-1) c(11) = 0
c†(11) c(1-1)     c†(11) c(2-1) + c†(11) c(3-1) - c†(21) c(1-1) - c†(31) c(1-1) = 0
c†(1-1) c(1-1)     c†(1-1) c(2-1) + c†(1-1) c(3-1) - c†(2-1) c(1-1) - c†(3-1) c(1-1) = 0
c†(1-1) c†(11) c(1-1)     4 c†(1-1) c†(11) c(1-1) + c†(1-1) c†(11) c(2-1) + c†(1-1) c†(11) c(3-1) - c†(1-1) c†(21) c(1-1) - c†(1-1) c†(31) c(1-1) + c†(11) c†(2-1) c(1-1) + c†(11) c†(3-1) c(1-1) = 0
c(1-1) c(11)     -4 c(1-1) c(11) + c(1-1) c(21) + c(1-1) c(31) - c(11) c(2-1) - c(11) c(3-1) = 0
c†(11) c(1-1) c(11)     -4 c†(11) c(1-1) c(11) + c†(11) c(1-1) c(21) + c†(11) c(1-1) c(31) - c†(11) c(11) c(2-1) - c†(11) c(11) c(3-1) - c†(21) c(1-1) c(11) - c†(31) c(1-1) c(11) = 0
c†(1-1) c(1-1) c(11)     -4 c†(1-1) c(1-1) c(11) + c†(1-1) c(1-1) c(21) + c†(1-1) c(1-1) c(31) - c†(1-1) c(11) c(2-1) - c†(1-1) c(11) c(3-1) - c†(2-1) c(1-1) c(11) - c†(3-1) c(1-1) c(11) = 0
c†(1-1) c†(11) c(1-1) c(11)     c†(1-1) c†(11) c(1-1) c(21) + c†(1-1) c†(11) c(1-1) c(31) - c†(1-1) c†(11) c(11) c(2-1) - c†(1-1) c†(11) c(11) c(3-1) - c†(1-1) c†(21) c(1-1) c(11) - c†(1-1) c†(31) c(1-1) c(11) + c†(11) c†(2-1) c(1-1) c(11) + c†(11) c†(3-1) c(1-1) c(11) = 0
```

可能出问题的地方包括
- 离谱的bug导致约束计算错误
- 老问题：虚部实部分离时乘法出错
- 太多变量没有出现在约束条件中？

观察是否各个约束都计算正确。以
```
c†(1-1) c†(11) c(1-1) c(11)  =>  c†(1-1) c†(11) c(1-1) c(21) + c†(1-1) c†(11) c(1-1) c(31) - c†(1-1) c†(11) c(11) c(2-1) - c†(1-1) c†(11) c(11) c(3-1) - c†(1-1) c†(21) c(1-1) c(11) - c†(1-1) c†(31) c(1-1) c(11) + c†(11) c†(2-1) c(1-1) c(11) + c†(11) c†(3-1) c(1-1) c(11) = 0
```
为例，在`2022-5-4.nb`中我们有
```
hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 1, "up"}, {c, 1, "down"}, {c, 2, "up"}] + 
 hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 1, "up"}, {c, 1, "down"}, {c, 3, "up"}] - 
 hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 1, "up"}, {c, 1, "up"}, {c, 2, "down"}] - 
 hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 1, "up"}, {c, 1, "up"}, {c, 3, "down"}] - 
 hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 2, "up"}, {c, 1, "down"}, {c, 1, "up"}] - 
 hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 3, "up"}, {c, 1, "down"}, {c, 1, "up"}] + 
 hubbardOpString[{SuperDagger[c], 1, "up"}, {SuperDagger[c], 2, "down"}, {c, 1, "down"}, {c, 1, "up"}] + 
 hubbardOpString[{SuperDagger[c], 1, "up"}, {SuperDagger[c], 3, "down"}, {c, 1, "down"}, {c, 1, "up"}]
```

| | |
| :------ | :------|
| `c†(1-1) c†(11) c(1-1) c(21)` | `hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 1, "up"}, {c, 1, "down"}, {c, 2, "up"}]`|
| ` c†(1-1) c†(11) c(1-1) c(31)` | `hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 1, "up"}, {c, 1, "down"}, {c, 3, "up"}]` |
| `- c†(1-1) c†(11) c(11) c(2-1)` | `-hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 1, "up"}, {c, 1, "up"}, {c, 2, "down"}]` |
| `- c†(1-1) c†(11) c(11) c(3-1)` | `- hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 1, "up"}, {c, 1, "up"}, {c, 3, "down"}]` |
| `- c†(1-1) c†(21) c(1-1) c(11)` | `- hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 2, "up"}, {c, 1, "down"}, {c, 1, "up"}]` |
| `- c†(1-1) c†(31) c(1-1) c(11)` | `- hubbardOpString[{SuperDagger[c], 1, "down"}, {SuperDagger[c], 3, "up"}, {c, 1, "down"}, {c, 1, "up"}]` |
| `c†(11) c†(2-1) c(1-1) c(11)` | `hubbardOpString[{SuperDagger[c], 1, "up"}, {SuperDagger[c], 2, "down"}, {c, 1, "down"}, {c, 1, "up"}]` |
| `c†(11) c†(3-1) c(1-1) c(11)` | `hubbardOpString[{SuperDagger[c], 1, "up"}, {SuperDagger[c], 3, "down"}, {c, 1, "down"}, {c, 1, "up"}]` |

应该是没问题了。

至少现在有一个问题，就是$K=5$时能够跑起来而$K=4$不行，说明$K=4$时，满足$l(O) \leq K / 2$的算符选两个乘起来可能会突破算符空间。这个是很奇怪的。

现在可以做的事情：
- 添加对称性的功能

关于粒子数算符和自旋算符的$\lang N, O \rang = 0$和$\lang S^z, O \rang = 0$约束，打印出来看一看：
```julia
N = sum(i -> cdag(i, 1) * c(i, 1) + cdag(i, -1) * c(i, -1), site_list) |> normal_form
Sz = sum(i -> cdag(i, 1) * c(i, 1) - cdag(i, -1) * c(i, -1), site_list) |> normal_form

for op in hubbard_opstr_basis
    println(comm(N, op) |> normal_form)
end

for op in hubbard_opstr_basis
    println(comm(Sz, op) |> normal_form)
end
```
发现输出从来不会多于一项——正确。

如下代码
```julia
count = 0
for op in hubbard_opstr_basis
    res = comm(N, op) |> normal_form
    if length(res.terms) == 0
        count += 1
    end
end
println(count)
```
给出110，`hubbard_opstr_basis`总共388个，然后`particle_number_constraint_ops`长度为278。一致。

我有点被整糊涂了，比如说site 2湮灭，然后site 1产生，这样的算符期望值是零吗？？？

DMRG的程序是容易写的，但是计算出来的能量和韩希之那篇Hubbard模型的文章中列举的参考数据不相符。在$t = 1, U = 4$时，我得到的单个site的能量是-0.81左右，但是韩希之的文章说是-0.58左右。
我猜测这是因为没有让总电子数保持为每个site分配一个电子的原因，把图画出来，是`2022-5-4.png`，可以看到明显每个site没有一个电子。
是否这是能量总体偏低的原因呢？

看来是的，需要加入能量守恒的条件。这么做了以后，我们得到`-56.999029759703596`的总能量，分到每个site上差不多-0.57，正好对上。

运行
```julia
plot(expect(ψ,"Nupdn"), legend = false)
```
得到`2022-5-4-2.png`。这和韩希之那篇文章上的结果也是一致的。

观察
```julia
 plot(correlation_matrix(ψ, "Cdagup", "Cup")[1, 1:20], legend = false)
```
会发现$c^\dag_i c_j$在$i \neq j$时也是可以有非零值的。

我们有
```
julia> expect(ψ, "Cdn")
100-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
 0.0
 ⋮
 0.0
 0.0
 0.0
 0.0
 0.0
```
因此看起来单独一个$c$算符的期望值确实是零。

```
julia> correlation_matrix(ψ, "Cdagdn", "Cdagdn")
100×100 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 ⋮                        ⋮                        ⋱  ⋮                        ⋮
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  
```
因此两个$c$乘积的期望值确实是零。

在`dmrg.jl`中频繁地出现一些错误：
- `ERROR: Parity-odd fermionic terms not yet supported by AutoMPO`
- `ERROR: No block found with QN equal to QN(("Nf",-2,-1),("Sz",0))`

```julia
opstr_mpo = OpSum()
opstr_mpo += "Cdagup", 50, "Cdagdn", 50, "Cdn", 50
inner(ψ, MPO(opstr_mpo, sites), ψ)
```
确实给出`ERROR: Parity-odd fermionic terms not yet supported by AutoMPO`。

```julia
opstr_mpo = OpSum()
opstr_mpo += "Cdagup", 50, "Cdagdn", 50, "Cdagdn", 50, "Cup", 49
inner(ψ, MPO(opstr_mpo, sites), ψ)
```
给出`No block found with QN equal to QN(("Nf",-2,-1),("Sz",2))`，虽然它显然是零！

```julia
opstr_mpo = OpSum()
opstr_mpo += "Cdagup", 50, "Cdagdn", 50, "Cdn", 50, "Cup", 49
inner(ψ, MPO(opstr_mpo, sites), ψ)
```
就不报错。

```julia
for op in particle_number_constraint_ops
    println(benchmark_point_dmrg[op])
end

##

for op in spin_constraint_ops
    println(benchmark_point_dmrg[op])
end
```
输出的结果全都是零，因此这些约束是正确的。

通过DMRG获得的期望值都是实数；在构造优化问题时，我们可以不必分开虚部和实部，直接上手做。原本区分虚部实部的代码放在`optimization_problem-original.jl`中，我们去处理`optimization_problem.jl`。

报奇怪的错，什么variable not owned之类的。

我们索性一次把benchmark做完吧。

```julia
(
    - t * benchmark_point_dmrg[cdag(1,-1) * c(2,-1)] 
    - t * benchmark_point_dmrg[cdag(1, -1) * c(3, -1)] 
    - t * benchmark_point_dmrg[cdag(1, 1) * c(2, 1)] 
    - t * benchmark_point_dmrg[cdag(1, 1) * c(3, 1)] 
    - U * benchmark_point_dmrg[cdag(1, -1) * cdag(1, 1) * c(1, -1) * c(1, 1)]
)
```
给出`-0.5736684912615906`，这个是正确的。

运行
```julia
hubbard_opstr_basis_expected_dmrg_values = map(hubbard_opstr_basis) do op
    benchmark_point_dmrg[op]
end

for constraint_coefficients in H_constraints_coefficients
    println(constraint_coefficients' * hubbard_opstr_basis_expected_dmrg_values)
end
```
得到
```
0.0     
0.0
0.0
0.0
-2.7755575615628914e-17
0.0
0.0
-4.718447854656915e-16
3.0531133177191805e-16
0.0
0.0
-4.996003610813204e-16
0.0
-3.608224830031759e-16
0.0
0.0
0.0
0.0
7.216449660063518e-16
0.0
0.0
0.0
0.0
```
因此最大的violation多半是数值舍入误差。

运行
```julia
M_dmrg = zeros(length(M_mat_spanning_opstr_indices), length(M_mat_spanning_opstr_indices))

for i in 1 : length(M_mat_spanning_opstr_indices)
    for j in 1 : length(M_mat_spanning_opstr_indices)
        M_dmrg[i, j] = M_coefficient[i, j]' * hubbard_opstr_basis_expected_dmrg_values
    end
end

eigen(M_dmrg).values
```
发现本征值好像都是正的。

这样benchmark应该是没问题的；但是，运行结果仍然是`Status: Dual_infeasible`。