复现2006.06002中对Hubbard模型的处理
======

Created on 2022.1.17.

我们这里尝试复现2006.06002中对Hubbard模型的处理，并且为更多的模型提供一个原型。

模型是
$$
H=-\sum_{\langle x y\rangle \sigma} c_{x \sigma}^{\dagger} c_{y \sigma}+U \sum_{x} n_{x \uparrow} n_{x \downarrow}
$$
在二维正方晶格上面。

# 理论

- Eq. (2)中的$C_\alpha$和$U_\alpha$包括：
  - 连续对称性：总粒子数、总自旋
  - 离散对称性：动量、旋转
- $\mathcal{C}_2$集合（保证正定性）
  - 由于局域性，由足够短的费米算符串张成
- $\mathcal{C}_1$集合（保证共轭、连续对称性生成元守恒、离散对称性成立）
  - 由$\mathcal{C}_2$中两个串相乘得到的全体算符张成
- 应当在Eq. (2)的约束下优化$F[H]$，这能够给系统总能量设置一个下界（是下界，因为真实的能量肯定在可行域中）

# 算法纲要

## $\mathcal{C}_2$的表示

为了避免double counting，应当始终使用正规序。

2022.1.17 各点的排序方式应该是怎样的？

## 泛函$\mathcal{F}$的表示

泛函$\mathcal{F}[O] := \langle \rho O \rangle$是一个线性泛函。按照前述$C_1$和$C_2$的构造，我们可以将$\mathcal{F}$看成线性映射$\mathcal{C}_2 \times \mathcal{C}_2 \to \Complex$。
最为朴素的实现是将$\mathcal{F}$当成一个列矢量，其基底是$\mathcal{C}_2 \times \mathcal{C}_2$。

Poor man's OPE: 用反对易关系把一个算符串重排成一系列正规序算符的线性组合然后代入$\mathcal{F}$中；一个费米算符出现两次就认为是零（因为不可能有两个一样的费米子）。

## 晶格

由于在无穷大的晶格上工作（远处的算符平移到原点附近来求值），不使用任何基于有限晶格的库。看起来，将每个格点都标上号是没有意义的，于是我们使用坐标$(x, y)$来标记格点。

哈密顿量中的各项也是平移到原点附近求值，最后整个哈密顿量正比于格点总数$N$，把$N$扔掉以后对能量做最优化。
$$
H=-\sum_{\langle x y\rangle \sigma} c_{x \sigma}^{\dagger} c_{y \sigma}+U \sum_{x} n_{x \uparrow} n_{x \downarrow}
$$
的第一项正比于$2N$，第二项正比于$N$，即
$$
\begin{aligned}
    \langle H \rangle &= - 2 N t \sum_\sigma \langle c^\dagger_{0 \sigma} c_{\hat{\boldsymbol{x}} \sigma} \rangle + U N \lang n_{0\uparrow} n_{0\downarrow} \rang  \\
    &= - 2 N t \sum_\sigma \langle c^\dagger_{0 \sigma} c_{\hat{\boldsymbol{x}} \sigma} \rangle + U N \lang c^\dagger_{0 \uparrow} c^\dagger_{0 \downarrow} c_{0 \downarrow} c_{0 \uparrow} \rang .
\end{aligned}
$$
类似的，
$$
\lang N \rangle = N \sum_\sigma \lang c^\dag_{0 \sigma} c_{0 \sigma} \rang,
$$
以及
$$
\lang S_z \rang = N (\lang c^\dag_{0 \uparrow} c_{0 \uparrow} \rang - \lang c^\dag_{0 \downarrow} c_{0 \downarrow} \rang).
$$

## 施加对称性

有待解决的问题：是否有把握
- 或者所有约束条件中涉及的算符序列都可以通过Poor man's OPE划归到$\mathcal{C}_2$中，或是
- 所有约束条件都能够化归为对$\mathcal{F}$的约束？看起来和$\mathcal{C}_1$有关的约束主要是用于构造$\mathcal{F}$的；主要的计算量都在正定性条件上面。

- 平移对称性：
- 旋转对称性：
- 复共轭：似乎可以将$\mathcal{C}_2$集合中互为共轭转置的算符配对

## 正定性

我怀疑这才是真的关键的地方……这里的约束看起来是不能够通过限制$\mathcal{F}$的形式引入的

2022.1.17 感觉应该先做一个toy model感觉一下情况；应该看一下2004.10212那篇文章

2022.1.30 目前的计划：
- 使用$\lang c^\dag_{i \sigma} c_{j \sigma} \rang$和$\lang n_{i \uparrow} n_{i \downarrow} \rang$做变量
- 优化目标$E$的表达式是显然的
- 通过$\lang H \mathcal{O} \rang$可以制造出适当的OPE
- 这样的问题是，无法在$\lang c^\dag_{i \sigma} c_{j \sigma} \rang$和$\lang n_{i \uparrow} n_{i \downarrow} \rang$之间建立联系……
  - 我需要一个Hubbard模型版本的位力定理，
  - 或者说要仅仅使用$\lang c^\dag_{i \sigma} c_{j \sigma} \rang$和$\lang n_{i \uparrow} n_{i \downarrow} \rang$构造一个方程

最好是拿mathematica算一算对易子

# 实验记录

## 2022.2.8

在`hubbard-simple-prototype\2022-2-8.nb`中计算费米子对易子。

下面要写一个程序弄一个类似于`hubbard-simple-prototype\2d-square-lattice-indexing.pdf`的格点编码。

## 2022.2.9

完成格点编码。准备构造$\mathcal{O}_2$集合。在`hubbard-simple-prototype\2022-2-9.nb`中做这件事。

## 2022.2.13

在`hubbard-simple-prototype\2022-2-13.nb`中继续构造$\mathcal{O}_2$集合。
撞上了（我预期中肯定会撞上的）问题，就是粗放地考虑的话，可能的算符太多了。Mathematica就没有迭代器什么的吗……

发现`comm`函数定义有问题；应该是修正了该问题。加入了`pauliRule`规则，用以去除$c^2$项。

~好像`comm`的定义还是有问题，就是计算结果不对易。~我脑子真的是进水了，有个🔨问题，没有问题。

## 2022.2.22

$\mathcal{C}_2$其实没有必要被显式地构造；反正唯一会用到$O \in \mathcal{C}_2$的地方就是$\mathcal{F}[O^\dagger O]$，保证是厄米的，那么索性只构造$O^\dagger O$的产生算符部分（这部分的程序我也已经写好了），湮灭算符部分靠dagger弄出来就行了；然后这样$l(O)$也可以自动计算（所以$l(O)$定义成这样大概是为了可加性）。

我可以先试试不去根据$l(O)$筛选算符，先看看最核心的bootstrap能不能做。

## 2022.3.29

戚老师的说法是，也许Hubbard model就是适合linear SDP，不需要特别高精度的库呢？

我感觉可以老实读一读那篇矩阵模型的文章了。它倒是没说什么特别的东西，不过无论如何算符序列都不很长。

目前在`operator-string.jl`中构造合格的算符串，虽然还没有跑通。

## 2022.4.2

~尝试使用最土鳖的冒泡排序法来解决算符normal ordering的问题。在`operators.jl`中做这件事。~
没有完成算符序列的构造没办法处理任意的算符啊。

连续两次运行`operator-string.jl`程序卡死，已经说明`qualified_opstr`是非常大的了。

发现`operator-string.jl`编写有误。发现以下错误：
- 计算算符串大小的时候没有把算符总数算进去
- 没有及时pop`current_op_str`

然后不知道为什么程序就一直运行着，哪怕我设置了最多迭代100次。尝试调试。然后撞上了老bug。

现在的问题是好像一直在1号和9号格点上一个up一个down的构型上打转……

现在`operator-string.jl`应该大体上没有问题了。需要进一步测试。

## 2022.4.3

把`operator-string.jl`中的`println`都注释掉了。

我有点不想检查是不是能够cover全部情况了……`operator-string.jl`中构造的算符只要是合法的，就能够参加bootstrap。
算符代数那些东西更加容不得差错。唉但是好像还是不对，因为还是需要使用向量来保存算符展开的系数，所以肯定还是要给每个算符串编上号。

## 2022.4.6

我感觉我好像搞错了什么。在之前关于$x$和$p$的问题中，$\lang O^\dagger \rangle = \lang O \rang^*$被考虑进去了没有？

## 2022.4.7

似乎有一种方法可以绕过（姑且这么说）在$\mathcal{C}_2$中构造非常复杂的算符空间的必要性。你想，(2)中的算符都是$O^\dag O$这种形式的，那么显然它们可以被normal order，然后其中的每一项的产生算符部分的$l$值当然都小于$K$，湮灭算符当然也是。
另一方面，我们有
> $\mathcal{C}_1$ spanned by the strings that appear in the products of two operators in $\mathcal{C}_2$

因此，看起来，可以使用类似于谐振子的方法，如下定义$M$矩阵和constraint：使用产生算符

这样就没有必要显式构造$\mathcal{C}_2$了，从而免去了将`operator-string.jl`的结果重新标上产生和湮灭标签的时间。
因此，只需要保证`operator-string.jl`的结果的完备性就可以了。或者其实完备性也用不着保证，反正如果计算过程中有算符落在$\mathcal{C}_1$外面就忽略这个计算。
好得很，好得很，，，

然后，既然现在某个`current_op_str`不再是一个$\mathcal{C}_1$算符中所有产生湮灭算符的标签，显然我们需要允许它完全没有经过原点。
平移对称性的事情怎么引入呢？不知道，但总之`operator-string.jl`不是讨论这个问题的地方。

现在的任务包括：
- 修改`operator-string.jl`使之能够容纳原点一次也没有出现的点序列
- 开发`operators.jl`

先做第一件事。将
```julia
current_op_str = [:up]
```
修改为
```julia
current_op_str = [:no]
```
看看会怎么样

## 2022.4.10

无论如何要开始好好实现了。

在`operator-string.jl`中运行
```julia
sort(values(hubbard_opstr_index) |> collect) == 1 : length(hubbard_opstr_basis)
```
得到`true`。

算费米子乘积的normal ordering似乎比我想象的要简单。

## 2022.4.11

今天需要做两件事：
- 一个是完全技术性的，就是让`operator-string.jl`中横行的全局变量包到`begin ... end`块里面
- 一个是使用QuantumAlgebra.jl来做费米子代数

第一件事似乎已经完成了；我感觉应该将`operator-string.jl`重命名为`operator-label.jl`，而将`operator.jl`重命名为`operator-algebra.jl`。
我也将`operator-string.jl`最后一个region挪到了`operator-algebra.jl`里面。

需要将`:up`替换成`:↑`。否则会报
```
ERROR: ArgumentError: index names must be single character, got up.
Stacktrace:
 [1] QuantumAlgebra.QuIndex(::String) at C:\Users\wujin\.julia\packages\QuantumAlgebra\6pQhc\src\operator_defs.jl:37
 [2] QuIndex at C:\Users\wujin\.julia\packages\QuantumAlgebra\6pQhc\src\operator_defs.jl:33 [inlined]     
 [3] iterate at .\generator.jl:47 [inlined]
 [4] collect_to!(::Array{QuantumAlgebra.QuIndex,1}, ::Base.Generator{Tuple{Int64,Symbol},Type{QuantumAlgebra.QuIndex}}, ::Int64, ::Int64) at .\array.jl:732
 [5] collect_to_with_first! at .\array.jl:710 [inlined]
 [6] collect(::Base.Generator{Tuple{Int64,Symbol},Type{QuantumAlgebra.QuIndex}}) at .\array.jl:691        
 [7] make_indices(::Int64, ::Symbol) at C:\Users\wujin\.julia\packages\QuantumAlgebra\6pQhc\src\operator_defs.jl:70
 [8] FermionDestroy(::QuantumAlgebra.QuOpName, ::Int64, ::Vararg{Any,N} where N) at C:\Users\wujin\.julia\packages\QuantumAlgebra\6pQhc\src\operator_defs.jl:119
 [9] c(::Int64, ::Vararg{Any,N} where N) at C:\Users\wujin\.julia\packages\QuantumAlgebra\6pQhc\src\operator_defs.jl:248
 [10] top-level scope at REPL[9]:1
```

记得要引人家的包！

注意到
```
julia> normal_form(c(1,:↑) * c(2,:↑))
c(1↑) c(2↑)

julia> normal_form(cdag(1,:↑) * cdag(2,:↑))
c†(1↑) c†(2↑)
```
如果我们要c†(1↑) c†(2↑) c(2↑) c(1↑)这种效果，可能需要去访问`QuExpr`的`items`成员然后手动把负号加上去。

目前，我们**不**采用这种做法，就是说直接自乘，不去保持c†(1↑) c†(2↑) c(2↑) c(1↑)。

现在算符代数部分写应该是写完了，对不对那就不知道了……
新增一个文件`configuration.jl`，用于保存各种参数。

`Construct the M matrix`一节是有问题的，这样只能够剩下
```
julia> M_mat_spanning_opstr_indices
1-element Array{Int64,1}:
 1
```

目前的问题是
```
KeyError: key c†(1-1) c†(11) not found
```

这个解决了。

然后处理等式约束。这个需要
- 根据$O$中涉及的算符构造哈密顿量
- 构造对称性（也许可以先放放？在这里和哈密顿量对易能够隐式施加对称性约束吗？）

## 2022.4.12

现在需要解决这个问题：计算算符和哈密顿量的对易关系时，hopping term要怎么构造？我感觉好像能够hop到的点都是在`site_list`里面。
就是你想，如果$l(O) \leq K$，那么一个格点的$1-$范数应该小于等于$K-1$，就是说其每个坐标的绝对值小于等于$K-1$，但是，`site_list`里面的点排成$(2K+1)^2$方阵。

`hopping-legs.pdf`展现了应该如何不double counting地写下hopping terms：找所有格点，写下每个格点向周围跑的hopping term，就可以（黄色方框标出了一个完整的hopping term $- t \sum_{\lang \boldsymbol{i}, \boldsymbol{j} \rang, \sigma} c^\dagger_{\boldsymbol{i} \sigma} c_{\boldsymbol{j} \sigma}$）

现在不知道为什么`H_hubbard`都是零。

## 2022.4.13

有待做的事情：
- 单个格点的哈密顿量
    $$
    \begin{aligned}
        H_1 &= -t \sum_{\sigma} c^\dagger_{1 \sigma} (c_{2 \sigma} + c_{4 \sigma} + c_{6 \sigma} + c_{8 \sigma}) + U n_{1 \uparrow} n_{1 \downarrow} 
    \end{aligned}
    $$

    注意这里没有h.c.项，即没有从1到2, 4, 6, 8号格点的跃迁，以免double counting。
- 平移对称性；实际上，需要先施加平移对称性，才能够保证以单格点哈密顿量为优化目标是合理的

## 2022.4.14

反之总之是写完了，错肯定是有的，先看看跑起来效果怎么样吧……

需要的文件：
- `main.jl`
- `configuration.jl`
- `operator-label.jl`
- `operator-algebra.jl`
- `optimization_problem.jl`
- `run_optimization.jl`
- `correlation-functions.jl`
- `2022-4-14-1.pbs`

核对
- PBS文件执行的确实是`main.jl`
- 命名正确

现在不知道为什么总是说`@match`找不到……

## 2022.4.15

现在可以做的事情包括：检测`@match`怎么回事。

## 2022.4.16

用`2022-4-16-run-1.sh`提交
- `main.jl`
- `configuration.jl`
- `operator-label.jl`
- `operator-algebra.jl`
- `optimization_problem.jl`
- `run_optimization.jl`
- `correlation-functions.jl`
- `2022-4-14-1.pbs`

依照报错修改代码，修改`optimization_problem.jl`。

不feasible可还行。

## 2022.4.20

所以benchmark的数据去哪里找呢……

为了benchmark，至少有一件事是可以做的，就是把`M`矩阵和约束都打印出来。