谐振子的数值bootstrap
======

考虑到谐振子的简单性我们还是从谐振子出发比较好……

# 理论

我们根据2004.10212中的Bootstraping the quantum anharmonic oscillator一节的做法来操作。
手工计算，从(2)出发，不断使用对易关系将左右两边所有的$x$都挪动到$p$前面，会得到
$$
\lang p^2 \rang = \lang x^2 \rang + 2 g \lang x^4 \rang, 
$$
于是就可以得到(3)。一个更加一般的版本是(4)。

注意$E$实际上就是$H$的期望值，我们有$\lang \mathcal{O} H \rang = E \lang \mathcal{O} \rang$。
这又导出了(5)。

从(4)和(5)出发可以推导出公式(6)。这个公式让我们能够递归地计算$\lang x^n \rang$。
于是，现在将用于bootstrap的算符集合选取为(7)中的集合（阶数从$0$到$K$的$x^i$的线性组合），即可得到一个形如$\bold{x}^\top \bold{M} \bold{x} \geq 0$的不等式，从而能够为$E$和$\lang x^2 \rang$确定范围。

# 算法纲要

我们有以下理论结果：
- $\lang x^0 \rang = 1$；
- $\lang x^2 \rang$是一个待优化的变量；
- $\lang x^4 \rang$由$E = 2 \lang x^2 \rang + 3 g \lang x^4 \rang$（原文公式(3)）确定，是一个待优化的变量（实际上我们是优化$E$）
- 高阶的$x^n$由公式(6)即
  $$
  4 t E\left\langle x^{t-1}\right\rangle +t(t-1)(t-2)\left\langle x^{t-3}\right\rangle -4(t+1)\left\langle x^{t+1}\right\rangle-4 g(t+2)\left\langle x^{t+3}\right\rangle=0
  $$
  给出；
- 奇数$x^n$期望值为零；
- 矩阵$M_{ij} = \lang x^{i+j} \rang$正定。

根据这些结果确定$E$和$\lang x^2 \rang$的范围。

2022.1.19: TODO 非简谐振子的微扰解，用来benchmark

# 细节

## 计算OPE的流程

- 先实现两个$x^m p^n$型算符的乘积的展开
- 对两个任意的输入算符$A$和$B$，遍历其稀疏矩阵中的非零元素，分别计算OPE以后加总

# 实验记录

## 解析计算

### 2022.1.19 

尝试使用MMA。内容保存于`2022-1-19.nb`中。大概能够画出正确的边界，但是存在两个问题：
- 首先阴影区域不正确
- 其次，与原文图1相比，左下方的弧线缺失，能量似乎没有下界

公式(7)中的$c_i$为虚部，似乎是左下方弧线缺失的一个可能的原因？

我擦我搞错了，矩阵正定不能拿行列式判断……

### 2022.1.20

在`2022-1-20.nb`中复现了图1第一张图和第二张图。严格的bootstrap确定能够将$E$确定得比较精确。

下面要做的是使用更加常规的方法求解非简谐振子。尝试薛定谔方程数值解，最朴素的那种格式。

在`operators.jl`中使用Δx=0.05算出来的能量本征值：

 0.04984325707680226        
 0.14921428524939048
 0.24795126811106247
 0.34604795716465764
 0.44349789463145856
 0.5402944016355128
 0.6364305654411193
 0.7318992256443275
 0.8266929592062982
 0.9208040642014521

使用Δx=0.1算出来的：

 0.09937101841211175
 0.2968386939280462
 0.4917317992198175
 0.6839968161582278
 0.8735763586427825
 1.0604086888571334
 1.2444271451800792
 1.4255594600017991
 1.6037269387993491
 1.7788434621950395

好吧我发现是我在$\laplacian$的定义中少除了因子$\Delta x^2$。现在对了。

发现：在取Δx为0.1, 0.05, 0.01时，计算得到$g=1$时的基态能量为
 1.390739671807073
 1.3919487648892392
 1.3923355279431426

前两个值小于bootstrap得到的1.3921989772554817，第三个值大于。

可能需要通过变分法来获得一个上界估计。

在`variation.jl`中却是发生了一件很奇怪的事情，使用更复杂的拟设，能量反而上升了。离谱。

## 2022.1.23

现在的问题是，如何能够更加“数值”地做bootstrap（在Hubbard模型里面估计就不能够用mma算了）。
使用2006.06002中的记号，在本例中，$\mathcal{C}_1$和$\mathcal{C}_2$一致，均选取为$x^n$。
泛函$\mathcal{F}[\cdot]$通过$E$和$\lang x^2 \rang$得到定义。
对称性条件仅有$\lang [\mathcal{O}, H] \rang = 0$这一个被使用了，并且我们使用了其强形式，即假定$\mathcal{F}$是某个能量本征态。

看起来，似乎可以尝试写一个最优化程序，从另一个角度重新计算本算例。然后问题就来了，就是我们怎么样在“矩阵正定”的约束下做最优化。

我们将在`convex-approach.jl`中做这件事。

目前的版本遇到的麻烦是，高阶$\lang x^n \rang$含有很多包含$E$和别的$\lang x^n \rang$乘在一起的项的约束。
似乎`Convex.jl`不支持这个东西。

Numerical Bootstrap in Quantum Mechanics 2108.11416似乎值得一读。

## 2022.1.31

在`operators.jl`中做对易子自动计算。而且好像还搞错了……

## 2022.2.1

算了我还是老老实实解析算吧……按照McCoy's formula（见[此处](https://en.wikipedia.org/wiki/Canonical_commutation_relation#Generalizations)）
$$
p^m x^n = x^n p^m + \sum_{k=1}^{\min(m, n)} \frac{(- \ii)^k n! m!}{k! (n-k)! (m-k)!} x^{n-k} p^{m-k}
$$
结果见`oscillator-simple-prototype\commutation.jl`。

## 2022.2.2

新建虚拟环境`optimization`，在其中安装了`https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/`。
单元测试出现问题。正在和开发者了解情况……

## 2022.2.6

单元测试问题解决。在`oscillator-simple-prototype\nonconvex-sdp-approach.jl`中尝试数值计算最优化。

眼下的问题似乎是`add_sd_constraint!`没有起到应有的效果，计算结束后`m_mat_minimizer`压根不是正定的。

为了解决此问题，在`nonconvex-sdp\my-task-oscillator-prototype.jl`中观察改变程序的哪一部分以后能够至少真的做SDP而不是简单的将能量降到零。
可能存在的bug：
- 不能正确处理递归？
  - 看起来存在除此以外的问题，因为将`K`取为2，此时在计算`expected_xn`时一次递归也没有使用，但是最后计算结果仍然没有保持正定性
- 要求保持正定的函数不能够调用自定义函数？
  - 改用
  ```julia
  sd_constraint((E, x²)) = [
    1.0       0        x² ;
    0         x²       0  ;
    x²        0        E - 2x²
  ]
  ```
  之后计算结果仍然没有保持正定性。
- 目标函数就是$E$，是不是需要耍一点什么花招……
  - 似乎使的，因为在将优化变量改为x²和x⁴之后，程序正确地识别出只需要将这两者都取为零即可；推测这是由于梯度下降法：目标函数就是其中一个优化变量的话梯度不会有任何变化

在获得了以上结论后，在`nonconvex-sdp\my-task-oscillator-prototype-2.jl`中重新尝试。
1. 尝试将优化变量改成x²和x⁴，其余不变——然后发现还是啥玩意都是零……
2. 为了分析是不是递归的原因，我们在`nonconvex-sdp\single-constraint-sdp-simple-test-1.jl`中测试一个不那么平凡，但是所有东西的定义都没有使用递归的例子，然后发现似乎正定性约束仍然没有起作用
3. 正在和作者联系……

## 2022.2.7

续上表。

1. 和作者联系并且更新了版本，现在`nonconvex-sdp\single-constraint-sdp-simple-test-1.jl`已经能够跑了，但是好像无法优化，就是说无法收敛。
2. 实际上并没有那么没法收敛，因为`nonconvex-sdp\2022-2-6.nb`中解析计算的结果是`x2 -> 0.666667, E0 -> 1.66667`，而另一方面，使用`nonconvex-sdp\single-constraint-sdp-simple-test-1.jl`没完全收敛的结果则是`[x4, x2]`=
   ```
    0.11111277779546458
    0.666669166649511
   ```
3. 我们尝试将迭代步数增大，来看是不是能够最终收敛。将
   ```julia
   options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=200))
   ```
   改成
   ```julia
   options = SDPBarrierOptions(sub_options=IpoptOptions(max_iter=400))
   ```
   仍然没有收敛……比较这里的结果1.6666766666854158和解析计算出来的1.6666667242418427，相对误差5.965465937786967e-6。还算能够接受。这里的问题似乎在于，Ipopt会让需要保持为正定的矩阵的最小本征值非常解决零但是不是零。我不知道这个是bug还是feature，就是说如果极小值点是优化边界确定的，那正常结果是收敛还是不是。

下一个需要解决的问题来自`oscillator-simple-prototype\nonconvex-sdp-approach.jl`。使用目前版本的约束矩阵，即
```julia
sd_constraint((E, x²)) = [expected_xn(E, x², i + j) for i in 0 : K, j in 0 : K]
```
会导致如下报错：
```
ERROR: LoadError: MethodError: no method matching getindex(::Nothing, ::Int64)
```
被标红的行为
```julia
result = optimize(model, alg, x0, options = options)
```
以及
```julia
sd_constraint((E, x²)) = [expected_xn(E, x², i + j) for i in 0 : K, j in 0 : K]
```
不是很确定这到底是怎么回事。数组推导是不能用还是？？

Hubbard模型的bootstrap还tm遥遥无期……理论推导先试试吧……

在`nonconvex-sdp\single-constraint-sdp-size-variable-mat-1.jl`中观察大小任意的约束矩阵是不是会出问题。
好像并没有什么问题，优化结果是完全正确的。那就奇了怪了，哪儿出问题了

在`nonconvex-sdp\single-constraint-sdp-size-variable-mat-2.jl`中尝试融合`nonconvex-sdp\single-constraint-sdp-size-variable-mat-1.jl`和`nonconvex-sdp\my-task-oscillator-prototype-2.jl`，看哪儿可能会出问题。

似乎已经发现了会出问题的地方：将
```julia
Mij(x⁴, x², n - 2) * x² + Mij(x⁴, x², n - 4)
```
替换成
```julia
t = n - 3
    # According to (6) in 2004.10212
    (4t * (2x² + 3 * g * x⁴) * Mij(x⁴, x², t - 1) + t * (t - 1) * (t - 2) * Mij(x⁴, x², t - 3) - 4 * (t + 1) * Mij(x⁴, x², t + 1)) / (4g * (t + 2))
```
就会有问题出现。在新的REPL中重新运行程序，同样的问题再次出现，说明不是来自不同优化任务的冲突。

尝试简化上面这段代码，看看什么feature被用到的时候会报错。

1. 替换它为
   ```julia
   t = n - 3
    # According to (6) in 2004.10212
    (4t * (2x² + 3 * g * x⁴) * Mij(x⁴, x², t - 1) + t * (t - 1) * (t - 2) * Mij(x⁴, x², t - 3) - 4 * (t + 1) * Mij(x⁴, x², t + 1)) 
   ```
  出错
2. 替换它为
   ```julia
   t = n - 4
    Mij(x⁴, x², t + 2) * x² + Mij(x⁴, x², t)
   ```
   不出错
3. 替换它为
   ```julia
   t = n - 4
    t * Mij(x⁴, x², t + 2) * x² + Mij(x⁴, x², t)
   ```
   出错
4. 替换它为
   ```julia
   (n - 4) * Mij(x⁴, x², n - 2) * x² + Mij(x⁴, x², n - 4)
   ```
   出错。错误定位。

## 2022.2.8

提了[issue](https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/issues/6)。
使用作者提供的暂时解决方案后，`nonconvex-sdp\single-constraint-sdp-size-variable-mat-3.jl`中没有出现错误；预期的结果出现了。

尝试在`nonconvex-sdp\my-task-oscillator-prototype-3.jl`中按照一样的方式解决这个问题。

可以跑起来但是没有得到预期的结果。奇了怪了。能量似乎比允许的更低，奇怪……

然后发现好像正定性约束又失效了。啊……

我除了我是个弱智以外真的是一无所知。为了算对易子我耗了两天在`oscillator-simple-prototype\commutation.jl`上面，然后一分钟不到的`oscillator-simple-prototype\2022-2-8.nb`就把问题完全解决了。
真的是弱智，反正也不需要把程序弄得多自动化，最优化算法八字没一撇，我是脑子进了什么水要拿julia写完整个程序？？？

## 2022.2.22

根据上次讨论的意见，观察分项添加约束能否绕过求解器的bug。在`nonconvex-sdp\my-task-oscillator-prototype-4.jl`中做这件事。

这成功地避开了[这个github issue](https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/issues/6)中的自动微分问题。
然而，semidefinite constraint再次失效了。

下面有几件事可以做：
- 更换求解器
  - 阅读韩希之的[代码](https://github.com/hanxzh94/matrix-bootstrap)
  - JuMP
- 寻找[这个github issue](https://github.com/JuliaNonconvex/NonconvexSemidefinite.jl/issues/6)中的自动微分漏洞
- 提`nonconvex-sdp\my-task-oscillator-prototype-4.jl`涉及的新的正定性失效的漏洞

尝试在`oscillator-simple-prototype\jump-oscillator-1.jl`中用JuMP做优化。
重命名`optimization`环境为`nonconvex`环境；新建`jump`环境用于存放JuMP有关内容。

装CDCS还tm没法预编译，，，改用CSDP

我敢肯定`oscillator-simple-prototype\jump-oscillator-1.jl`中有东西写错了（约束位置不对，等等），不过我先跑一遍看看吧。

可以发现，CSDP支持linear SDP，不支持nonlinear constraints；反之，Ipopt支持nonlinear constraints，但是不支持SDP。

看起来，使用JuMP可能能写出更加清楚的代码（并且更不容易因为自动微分的问题而惹上麻烦），但是靠直接调包是做不成这个问题的。
阅读韩希之的代码似乎无法避免，，，

## 2022.2.24

应该好好读读Supplementary Material for 'Bootstrapping Matrix Quantum Mechanics' in 2004.10212.

2006.06002声称（"can be solved efficiently and accurately by SDP e.g. with [19, 20]"）

看起来，作者在做数值优化的时候只用到了约束$\lang [H, \mathcal{O}] \rang = 0$，因此只需要处理线性约束？

我们尝试写下一个基于上述做法的程序，在`oscillator-simple-prototype\jump-oscillator-2.jl`中。代码肯定是需要大大重写。

## 2022.2.25

现在的问题是怎么在同时考虑$x$和$p$的情况下构造$\mathcal{M}$矩阵；这里的问题是，generally speaking这是一个复矩阵，然后不知道CSDP支不支持就是了，，，

## 2022.3.2

安装`OffsetArrays.jl`。我实在是被这些什么0开头还是1开头的东西搞麻了。

我之前在`jump-oscillator-2.jl`中写的东西看起来是有很大问题的，就是我没有分清楚$M_{ij}$的指标对应的算符和出现在矩阵元上的算符的区别。

设$K$为$i$指标能够取的最大值，无论它表示什么，那么$M_{KK} = \lang O^\dag_K O_K \rang$是参与bootstrap的最高阶的算符了。对算符的编码必须编码到$O^\dag_K O_K$。

所以需要怎么设计编码呢？

## 2022.3.3

尝试使用`OffsetArray`来解决麻烦的从0还是从1开始的问题。

照惯例，又出问题了。为了保证计算对易子时不会出现index out of bound的问题，实际上需要将`L_max`取为参与bootstrap的算符长度加4（因为哈密顿量中最高出现了$x^4$项）；然后另外一方面，在协方差矩阵$M$中出现的$O^\dag O$形式的算符如果正好覆盖所有参与bootstrap而具有正定性的算符，那么$O$（从而，协方差矩阵的指标对应的算符）的长度就大体上应该是`(L_max - 4) / 2`。

所以这里头实际上有两个所谓的“index”：其一是算符系数向量的index，就是现在已经在`oscillator-simple-prototype\jump-oscillator-2.jl`里面定义的index，然后我们还需要M矩阵的那种index。

一时半会估计写不完，先随便弄点toy model吧。见`jump-toy-1.jl`

在`jump-toy-1-benchmark.nb`和`jump-toy-2-benchmark.nb`中玩了几把，好像SDP加线性约束还就是凸优化。
[这篇博客文章](https://blog.csdn.net/sdgyfbtnyj/article/details/100101233)反正是这么说的；到时候可以读一读优化方面的书。

[维基](https://en.wikipedia.org/wiki/Semidefinite_programming)本来应该好好读一读的。

`jump-toy-1.jl`和`jump-toy-2.jl`均测试成功。

## 2022.3.5

今天的讨论结果：
- 应当注意，linear SDP方法不能够给出清晰的第一激发态！（linear SDP给出的是凸的可行域，凸的可行域不可能分成一块一块）
- 是否可以通过$\lang O^\dagger O \rang \geq 0$导出不确定性关系？是否所有仅仅由算符代数确定的自洽性条件都可以使用某个算符的正定性条件给出？
  - 是否正定性条件普遍成立$\Leftrightarrow$我们的$\lang \cdot \rang$对应正定的密度矩阵？？

尝试在Julia里面找一个包做x和p的对易关系。QuantumAlgebra算Hubbard模型时肯定会用上，但是它好像不支持x和p。

在`2022-3-5.nb`中构造$M$矩阵，用linear SDP。从M[[27, 16]]的计算结果可以看出，最终构造出来的$M$矩阵是含有大量$\ii$的。
因此，我们需要确保优化器能够处理复变量问题。

在`jump-toy-3.jl`和`2022-3-5.nb`中讨论该问题。`2022-3-5.nb`表明Mathematica可以处理这类问题。`jump-toy-3.jl`表明JuMP不能直接处理这类问题。
事实上，我们有
```
julia> im * x1
ERROR: InexactError: Float64(im)
Stacktrace:
 [1] Real at .\complex.jl:37 [inlined]
 [2] convert at .\number.jl:7 [inlined]
```
一个可能的方法是将$\ii$替换成`[0 -1; 1 0]`。我们在`jump-toy-4.jl`中做这件事。成功了，和`2022-3-5.nb`完全一致。将有关内容放在`jump-toy-3-benchmark.nb`中。

## 2022.3.9

`jump-oscillator-2.jl`是写完了，但是报奇怪的错误。
```
LoadError: ArgumentError: Empty constraint MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.EqualTo{Float64}}(1): MathOptInterface.ScalarAffineFunction{Float64}(MathOptInterface.ScalarAffineTerm{Float64}[], 0.0)-in-MathOptInterface.EqualTo{Float64}(0.0). Not supported by CSDP.
```

## 2022.3.12

将`L_max`改成4，然后错误当然还是来了。这个错误真的很奇怪。要怎么排查呢？
- 首先至少应该看出来被报错的这个所谓empty constraint到底是什么
- 不知道从错误信息里面能不能看出来什么东西。
- 也许现在可以把约束都打印出来看一看

```
Empty constraint MathOptInterface.ConstraintIndex{
    MathOptInterface.ScalarAffineFunction{Float64},
    MathOptInterface.EqualTo{Float64}
}(1): MathOptInterface.ScalarAffineFunction{Float64}(MathOptInterface.ScalarAffineTerm{Float64}[], 0.0)-in-MathOptInterface.EqualTo{Float64}(0.0)
```
感觉这个`MathOptInterface.ScalarAffineTerm{Float64}[], 0.0`很有些问题。哪个施加约束的地方会有`[]`出现呢？
首先肯定不是SDP那个约束。

执行如下命令：
```julia
x_power=1; p_power=1
op = xpopstr_xp_power(x_power, p_power)
cons = comm_with_ham(op)
lhs = transpose(real(cons)) * xpopstr_basis_real + transpose(imag(cons)) * xpopstr_basis_imag 
```
发现`lhs`的`[1, 1]`和`[2, 2]`分量都是零，所以会有一个`0 == 0`的约束。是否这是问题呢？

有趣，`jump-toy-4.jl`中还确实没有出现过这种类型的约束。hmm……

为了判断是不是这样，在`jump-toy-4-2022-3-12-test-1.jl`中测试加入一个`0 == 0`的约束会怎么样。

好家伙果然报error了。看来CSDP不怕重复约束，但是怕空约束。

修改了
```julia
# Building constraints
# make sure the size of OH does not cause any out of bound error.
# Then imposing constraints to the variables
```
后面的那一段，但是好像还是有同样的error。

## 2022.3.13

在执行`jump-oscillator-2.jl`之后，执行如下代码
```julia
for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = transpose(real(cons))
    cons_imag = transpose(imag(cons))
    lhs = cons_real * xpopstr_basis_real + cons_imag * xpopstr_basis_imag
    if cons_real == cons_imag == zero_xpopstr
        continue
    end
    if cons_real == zero_xpopstr
        #@constraint(model, lhs[1, 2] == 0.0)
        println(x_power, " ", p_power)
        println(lhs[1, 2])
        println()
    elseif cons_imag == zero_xpopstr
        #@constraint(model, lhs[1, 1] == 0.0)
        println(x_power, " ", p_power)
        println(lhs[1, 1])
        println()
    else
        #@constraint(model, lhs .== O22)
        println(x_power, " ", p_power)
        println(lhs)
        println()
    end
end
```
发现的确有一些LHS的矩阵存在全零分量。这是怎么回事？进一步，执行
```julia
for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = transpose(real(cons))
    cons_imag = transpose(imag(cons))
    lhs = cons_real * xpopstr_basis_real + cons_imag * xpopstr_basis_imag
    if cons_real == cons_imag == zero_xpopstr
        continue
    end
    if cons_real == zero_xpopstr
        #@constraint(model, lhs[1, 2] == 0.0)
        println(x_power, " ", p_power)
        println("real part == 0")
        println(lhs[1, 2])
        println()
    elseif cons_imag == zero_xpopstr
        #@constraint(model, lhs[1, 1] == 0.0)
        println(x_power, " ", p_power)
        println("imag part == 0")
        println(lhs[1, 1])
        println()
    else
        #@constraint(model, lhs .== O22)
        println(x_power, " ", p_power)
        println(lhs)
        println()
    end
end
```
发现前三个选择分支从来没有被进入过。见鬼，`cons_real`和`cons_imag`是转置过的。

做出对应修改之后程序至少是能够跑起来了。

在`L_max = 5`时能量优化出来是`0.5921578324750953`。当然这个肯定很不准。

结果把`L_max`改成12又跑不了了。改成6都跑不了。

重复前述诊断，运行
```julia
for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = transpose(real(cons))
    cons_imag = transpose(imag(cons))
    lhs = cons_real * xpopstr_basis_real + cons_imag * xpopstr_basis_imag
    if cons_real == cons_imag == zero_xpopstr'
        continue
    end
    if cons_real == zero_xpopstr'
        #@constraint(model, lhs[1, 2] == 0.0)
        println(x_power, " ", p_power)
        println("real part == 0")
        println(lhs[1, 2])
        println()
    elseif cons_imag == zero_xpopstr'
        #@constraint(model, lhs[1, 1] == 0.0)
        println(x_power, " ", p_power)
        println("imag part == 0")
        println(lhs[1, 1])
        println()
    else
        #@constraint(model, lhs .== O22)
        println(x_power, " ", p_power)
        println(lhs)
        println()
    end
end
```
似乎无法定位出`0 == 0`约束。实际上这回的错误信息是
```
ERROR: LoadError: ArgumentError: Empty constraint MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.EqualTo{Float64}}(28): MathOptInterface.ScalarAffineFunction{Float64}(MathOptInterface.ScalarAffineTerm{Float64}[], 0.0)-in-MathOptInterface.EqualTo{Float64}(-2.0). Not supported by CSDP.
```
这个`-2.0`就很神奇，我什么时候加入过这种约束了。

最坏的情况当然就是有约束彼此冲突了，就是说我之前自动生成约束的程序有问题……

首先定位这个出问题的约束到底是哪里引进的。
```julia
@constraint(model, M[2i - 1 : 2i, 2j - 1 : 2j] .== real_part + imag_part)
```
是这个约束引入的吗？注释掉这个约束，发现错误仍然存在。不注释这个约束，而注释前面$\lang [O, H] \rang = 0$的约束没有错误。因此问题不在这里。

一种检查引入的约束有无互相冲突的办法是计算行列式。

想到一个稍微容易一些的办法。建立一个具有一模一样的变量的模型，然后一条一条往上面加约束，看看会怎么样。
目前`L_max = 6`，然后运行
```julia
for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = transpose(real(cons))
    cons_imag = transpose(imag(cons))
    lhs = cons_real * xpopstr_basis_real + cons_imag * xpopstr_basis_imag
    if cons_real == cons_imag == zero_xpopstr'
        continue
    end
    if cons_real == zero_xpopstr'
        #@constraint(model, lhs[1, 2] == 0.0)
        println(x_power, " ", p_power)
        println("real part == 0")
        println(lhs[1, 2])
        println()
    elseif cons_imag == zero_xpopstr'
        #@constraint(model, lhs[1, 1] == 0.0)
        println(x_power, " ", p_power)
        println("imag part == 0")
        println(lhs[1, 1])
        println()
    else
        #@constraint(model, lhs .== O22)
        println(x_power, " ", p_power)
        println(lhs[1, 1])
        println(lhs[1, 2])
        println(lhs[2, 1])
        println(lhs[2, 2])
        println()
    end
end
```
输出为
```
0 1     
real part == 0
2 xpopstr_expected[26] + 4 xpopstr_expected[78]

0 2
-12 xpopstr_expected[51] - 2
4 xpopstr_expected[28] + 8 xpopstr_expected[80]
-4 xpopstr_expected[28] - 8 xpopstr_expected[80]
-12 xpopstr_expected[51] - 2

0 3
-6 xpopstr_expected[1] - 36 xpopstr_expected[53]
-24 xpopstr_expected[26] + 6 xpopstr_expected[30] + 12 xpopstr_expected[82]
24 xpopstr_expected[26] - 6 xpopstr_expected[30] - 12 xpopstr_expected[82]
-6 xpopstr_expected[1] - 36 xpopstr_expected[53]

0 4
-12 xpopstr_expected[3] - 72 xpopstr_expected[55] + 24
-96 xpopstr_expected[28] + 8 xpopstr_expected[32] + 16 xpopstr_expected[84]
96 xpopstr_expected[28] - 8 xpopstr_expected[32] - 16 xpopstr_expected[84]
-12 xpopstr_expected[3] - 72 xpopstr_expected[55] + 24

1 0
real part == 0
-2 xpopstr_expected[2]

1 1
real part == 0
-2 xpopstr_expected[4] + 2 xpopstr_expected[52] + 4 xpopstr_expected[104]

1 2
-2 xpopstr_expected[25] - 12 xpopstr_expected[77]
-2 xpopstr_expected[6] + 4 xpopstr_expected[54] + 8 xpopstr_expected[106]
2 xpopstr_expected[6] - 4 xpopstr_expected[54] - 8 xpopstr_expected[106]
-2 xpopstr_expected[25] - 12 xpopstr_expected[77]

1 3
-6 xpopstr_expected[27] - 36 xpopstr_expected[79]
-2 xpopstr_expected[8] - 24 xpopstr_expected[52] + 6 xpopstr_expected[56] + 12 xpopstr_expected[108]       
2 xpopstr_expected[8] + 24 xpopstr_expected[52] - 6 xpopstr_expected[56] - 12 xpopstr_expected[108]        
-6 xpopstr_expected[27] - 36 xpopstr_expected[79]

1 4
24 xpopstr_expected[25] - 12 xpopstr_expected[29] - 72 xpopstr_expected[81]
-2 xpopstr_expected[10] - 96 xpopstr_expected[54] + 8 xpopstr_expected[58] + 16 xpopstr_expected[110]      
2 xpopstr_expected[10] + 96 xpopstr_expected[54] - 8 xpopstr_expected[58] - 16 xpopstr_expected[110]       
24 xpopstr_expected[25] - 12 xpopstr_expected[29] - 72 xpopstr_expected[81]

2 0
2
-4 xpopstr_expected[28]
4 xpopstr_expected[28]
2

2 1
2 xpopstr_expected[1]
-4 xpopstr_expected[30] + 2 xpopstr_expected[78] + 4 xpopstr_expected[130]
4 xpopstr_expected[30] - 2 xpopstr_expected[78] - 4 xpopstr_expected[130]
2 xpopstr_expected[1]

2 2
2 xpopstr_expected[3] - 2 xpopstr_expected[51] - 12 xpopstr_expected[103]
-4 xpopstr_expected[32] + 4 xpopstr_expected[80] + 8 xpopstr_expected[132]
4 xpopstr_expected[32] - 4 xpopstr_expected[80] - 8 xpopstr_expected[132]
2 xpopstr_expected[3] - 2 xpopstr_expected[51] - 12 xpopstr_expected[103]

2 3
2 xpopstr_expected[5] - 6 xpopstr_expected[53] - 36 xpopstr_expected[105]
-4 xpopstr_expected[34] - 24 xpopstr_expected[78] + 6 xpopstr_expected[82] + 12 xpopstr_expected[134]      
4 xpopstr_expected[34] + 24 xpopstr_expected[78] - 6 xpopstr_expected[82] - 12 xpopstr_expected[134]       
2 xpopstr_expected[5] - 6 xpopstr_expected[53] - 36 xpopstr_expected[105]

2 4
2 xpopstr_expected[7] + 24 xpopstr_expected[51] - 12 xpopstr_expected[55] - 72 xpopstr_expected[107]       
-4 xpopstr_expected[36] - 96 xpopstr_expected[80] + 8 xpopstr_expected[84] + 16 xpopstr_expected[136]      
4 xpopstr_expected[36] + 96 xpopstr_expected[80] - 8 xpopstr_expected[84] - 16 xpopstr_expected[136]       
2 xpopstr_expected[7] + 24 xpopstr_expected[51] - 12 xpopstr_expected[55] - 72 xpopstr_expected[107]       
```

这其中
```
2 0
2
-4 xpopstr_expected[28]
4 xpopstr_expected[28]
2
```
这一节比较神奇啊。第一个应该排除的问题就是，以下这段代码
```julia
for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = transpose(real(cons))
    cons_imag = transpose(imag(cons))
    lhs = cons_real * xpopstr_basis_real + cons_imag * xpopstr_basis_imag
    if cons_real == cons_imag == zero_xpopstr'
        continue
    end
    if cons_real == zero_xpopstr'
        #@constraint(model, lhs[1, 2] == 0.0)
    elseif cons_imag == zero_xpopstr'
        #@constraint(model, lhs[1, 1] == 0.0)
    else
        #@constraint(model, lhs .== O22)
    end
end
```
是错误的：应该让`cons`整体乘以`xpopstr_basis_real + xpopstr_basis_imag`，否则会漏掉诸如$(1 + \ii) (\Re \lang x p^3 \rang + \ii \Im \lang x p^3 \rang)$的式子的交叉项。
按照这个方法做了以后在`L_max = 6`时程序能够跑起来了，但是不收敛。输出如下：
```
CSDP 6.2.0
Iter:  0 Ap: 0.00e+000 Pobj:  0.0000000e+000 Ad: 0.00e+000 Dobj:  0.0000000e+000
Iter:  1 Ap: 8.94e-001 Pobj:  2.9135765e-003 Ad: 7.96e-001 Dobj: -1.6691110e+003
Iter:  2 Ap: 7.23e-001 Pobj: -2.1461302e-002 Ad: 7.47e-001 Dobj: -4.3513612e+003
Iter:  3 Ap: 5.69e-001 Pobj: -2.6843704e-002 Ad: 6.25e-001 Dobj: -4.5970090e+003
Iter:  4 Ap: 7.18e-001 Pobj: -3.3056704e-002 Ad: 6.05e-001 Dobj: -3.5877049e+003
Iter:  5 Ap: 5.88e-001 Pobj: -3.3587080e-002 Ad: 4.74e-001 Dobj: -3.6026251e+003
Iter:  6 Ap: 5.45e-001 Pobj: -3.9891326e-002 Ad: 4.98e-001 Dobj: -6.8414694e+003 
Iter:  7 Ap: 3.25e-001 Pobj: -1.9793391e-003 Ad: 4.17e-001 Dobj: -1.7296038e+004 
Iter:  8 Ap: 3.56e-001 Pobj: -1.4929152e-002 Ad: 2.59e-001 Dobj: -3.0082726e+004 
Iter:  9 Ap: 2.07e-001 Pobj: -3.9899597e-003 Ad: 1.69e-001 Dobj: -5.1484112e+004 
Iter: 10 Ap: 1.75e-001 Pobj: -6.8952358e-004 Ad: 2.04e-001 Dobj: -1.0369714e+005 
Iter: 11 Ap: 4.30e-002 Pobj: -9.3357851e-003 Ad: 2.62e-002 Dobj: -1.5779293e+005 
Iter: 12 Ap: 7.45e-002 Pobj: -2.7797489e-002 Ad: 3.78e-002 Dobj: -2.2335467e+005 
Iter: 13 Ap: 8.01e-002 Pobj: -2.5124854e-002 Ad: 5.10e-002 Dobj: -3.1249471e+005 
Iter: 14 Ap: 3.09e-002 Pobj: -2.0653000e-002 Ad: 2.19e-002 Dobj: -3.9387359e+005 
Iter: 15 Ap: 5.86e-002 Pobj: -3.0493792e-002 Ad: 6.99e-002 Dobj: -7.5747103e+005 
Iter: 16 Ap: 2.07e-002 Pobj: -4.0745847e-002 Ad: 1.48e-002 Dobj: -8.6652382e+005 
Iter: 17 Ap: 6.72e-002 Pobj: -6.7351231e-002 Ad: 3.49e-002 Dobj: -1.1878179e+006 
Iter: 18 Ap: 2.99e-002 Pobj: -5.5795053e-002 Ad: 2.49e-002 Dobj: -1.5524760e+006 
Iter: 19 Ap: 5.52e-002 Pobj: -4.8554588e-002 Ad: 2.14e-002 Dobj: -1.8574070e+006 
Iter: 20 Ap: 6.89e-002 Pobj: -9.7445510e-003 Ad: 4.81e-002 Dobj: -2.3399780e+006 
Iter: 21 Ap: 7.61e-002 Pobj:  4.4322528e-003 Ad: 7.10e-002 Dobj: -2.9793394e+006 
Iter: 22 Ap: 3.80e-002 Pobj:  2.1248791e-002 Ad: 5.50e-002 Dobj: -3.9235038e+006 
Iter: 23 Ap: 1.77e-001 Pobj:  6.1237672e-002 Ad: 1.05e-001 Dobj: -4.8925425e+006 
Iter: 24 Ap: 2.25e-001 Pobj:  6.6629249e-002 Ad: 1.86e-001 Dobj: -4.8298094e+006 
Iter: 25 Ap: 2.95e-001 Pobj:  6.6920334e-002 Ad: 3.25e-001 Dobj: -4.6516071e+006 
Iter: 26 Ap: 4.06e-001 Pobj:  8.2055383e-002 Ad: 3.41e-001 Dobj: -4.2925248e+006 
Iter: 27 Ap: 2.57e-001 Pobj:  1.3668061e-001 Ad: 2.40e-001 Dobj: -3.9866845e+006 
Iter: 28 Ap: 2.00e-001 Pobj:  2.4420053e-001 Ad: 1.82e-001 Dobj: -3.7464056e+006 
Iter: 29 Ap: 1.39e-001 Pobj:  4.3827913e-001 Ad: 1.29e-001 Dobj: -3.5838074e+006 
Iter: 30 Ap: 8.74e-002 Pobj:  7.6671766e-001 Ad: 9.60e-002 Dobj: -3.4702762e+006 
Iter: 31 Ap: 5.80e-002 Pobj:  1.3437788e+000 Ad: 7.02e-002 Dobj: -3.3895233e+006 
Iter: 32 Ap: 3.80e-002 Pobj:  2.2832417e+000 Ad: 5.53e-002 Dobj: -3.3249359e+006 
Iter: 33 Ap: 2.62e-002 Pobj:  3.7491887e+000 Ad: 4.71e-002 Dobj: -3.2670520e+006 
Iter: 34 Ap: 2.25e-002 Pobj:  6.4019780e+000 Ad: 3.88e-002 Dobj: -3.2168456e+006 
Iter: 35 Ap: 2.14e-002 Pobj:  1.1798876e+001 Ad: 2.60e-002 Dobj: -3.1839252e+006 
Iter: 36 Ap: 1.63e-002 Pobj:  2.0920682e+001 Ad: 2.41e-002 Dobj: -3.1553522e+006 
Iter: 37 Ap: 2.25e-002 Pobj:  4.7869841e+001 Ad: 2.48e-002 Dobj: -3.1267838e+006 
Iter: 38 Ap: 2.55e-003 Pobj:  5.6466862e+001 Ad: 2.19e-002 Dobj: -3.1026234e+006 
Iter: 39 Ap: 9.74e-004 Pobj:  6.1266697e+001 Ad: 2.17e-002 Dobj: -3.0829492e+006 
Iter: 40 Ap: 1.33e-003 Pobj:  7.0662482e+001 Ad: 1.62e-002 Dobj: -3.0715400e+006 
Iter: 41 Ap: 5.09e-004 Pobj:  7.5778772e+001 Ad: 1.92e-002 Dobj: -3.0520063e+006 
Iter: 42 Ap: 1.06e-003 Pobj:  8.9753980e+001 Ad: 1.23e-002 Dobj: -3.0333867e+006 
Iter: 43 Ap: 5.06e-004 Pobj:  9.8387121e+001 Ad: 1.43e-002 Dobj: -3.0107974e+006 
Iter: 44 Ap: 4.97e-004 Pobj:  1.0834594e+002 Ad: 1.98e-002 Dobj: -2.9840187e+006 
Iter: 45 Ap: 1.14e-003 Pobj:  1.3834944e+002 Ad: 9.39e-003 Dobj: -2.9712425e+006 
Iter: 46 Ap: 1.90e-004 Pobj:  1.4564973e+002 Ad: 1.09e-002 Dobj: -2.9516628e+006 
Iter: 47 Ap: 5.05e-004 Pobj:  1.6620453e+002 Ad: 1.00e-002 Dobj: -2.9325857e+006 
Iter: 48 Ap: 6.35e-004 Pobj:  1.9498984e+002 Ad: 1.15e-002 Dobj: -2.9145314e+006 
Iter: 49 Ap: 1.12e-003 Pobj:  2.6226669e+002 Ad: 1.39e-002 Dobj: -2.9014456e+006 
Iter: 50 Ap: 3.66e-004 Pobj:  3.0258776e+002 Ad: 5.65e-003 Dobj: -2.8942541e+006 
Iter: 51 Ap: 1.62e-004 Pobj:  3.2584136e+002 Ad: 4.04e-003 Dobj: -2.8875372e+006 
Iter: 52 Ap: 2.19e-004 Pobj:  3.6017357e+002 Ad: 9.64e-003 Dobj: -2.8700087e+006 
Iter: 53 Ap: 2.58e-004 Pobj:  4.0442632e+002 Ad: 4.00e-003 Dobj: -2.8639588e+006 
Iter: 54 Ap: 2.45e-004 Pobj:  4.5327880e+002 Ad: 5.48e-003 Dobj: -2.8561048e+006 
Iter: 55 Ap: 2.88e-004 Pobj:  5.2212742e+002 Ad: 5.76e-003 Dobj: -2.8469554e+006 
Iter: 56 Ap: 1.67e-004 Pobj:  5.7095875e+002 Ad: 5.03e-003 Dobj: -2.8384368e+006 
Iter: 57 Ap: 3.19e-004 Pobj:  6.7872678e+002 Ad: 6.73e-003 Dobj: -2.8293572e+006
Iter: 58 Ap: 1.35e-004 Pobj:  7.3926868e+002 Ad: 3.11e-003 Dobj: -2.8250135e+006
Iter: 59 Ap: 2.90e-004 Pobj:  8.9170553e+002 Ad: 5.14e-003 Dobj: -2.8212096e+006
Iter: 60 Ap: 7.53e-005 Pobj:  9.4707469e+002 Ad: 4.71e-003 Dobj: -2.8168009e+006
Iter: 61 Ap: 2.20e-004 Pobj:  1.1486775e+003 Ad: 5.63e-003 Dobj: -2.8138429e+006
Lack of progress.  Giving up!
Failure: return code is 7
Primal objective value: 2.2832417e+000
Dual objective value: -3.3249359e+006
Relative primal infeasibility: 5.61e-002
Relative dual infeasibility: 3.86e+000
Real Relative Gap: -1.00e+000
XZ Relative Gap: 3.94e+000
DIMACS error measures: 8.05e-002 0.00e+000 6.66e+000 0.00e+000 -1.00e+000 3.94e+000
objective_value(model) = -2.283241739287405
-2.283241739287405
```

不收敛问题在`L_max = 4`时也存在。

好像存在另一个问题，就是
```julia
@objective(model, Min, xpopstr_expected[xpopstr_index(2, 0)] + xpopstr_expected[xpopstr_index(0, 2)] + g * xpopstr_expected[xpopstr_index(4, 0)])
```
仍然是原来虚、实部混在一起的版本。换成
```julia
@objective(model, Min, 
    xpopstr_expected[xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :real)] + 
    xpopstr_expected[xpopstr_expected_real_imag_parts(xpopstr_index(0, 2), :real)] + 
    g * xpopstr_expected[xpopstr_expected_real_imag_parts(xpopstr_index(4, 0), :real)])
```
不工作，有invalid index error。卧槽原来`xpopstr_expected_real_imag_parts`直接输出变量而不是index。
去掉`xpopstr_expected[]`以后程序可以跑起来，但是memory allocation failed。（此时`L_max = 12`）

现在将`L_max`调小到6，然后还是不收敛，不过现在至少objective value是68.06036356465847。输出如下：
```
CSDP 6.2.0
Iter:  0 Ap: 0.00e+000 Pobj:  0.0000000e+000 Ad: 0.00e+000 Dobj:  0.0000000e+000 
Iter:  1 Ap: 8.94e-001 Pobj:  1.0393856e-001 Ad: 7.96e-001 Dobj: -1.6690661e+003 
Iter:  2 Ap: 7.23e-001 Pobj:  6.1576793e-001 Ad: 7.47e-001 Dobj: -4.3514769e+003 
Iter:  3 Ap: 5.69e-001 Pobj:  3.9555074e-001 Ad: 6.25e-001 Dobj: -4.5971255e+003 
Iter:  4 Ap: 7.18e-001 Pobj:  3.2698353e-001 Ad: 6.05e-001 Dobj: -3.5877331e+003 
Iter:  5 Ap: 5.88e-001 Pobj:  3.8394091e-001 Ad: 4.74e-001 Dobj: -3.6026388e+003 
Iter:  6 Ap: 5.46e-001 Pobj:  2.8370552e-001 Ad: 4.98e-001 Dobj: -6.8415674e+003 
Iter:  7 Ap: 3.24e-001 Pobj:  1.7443035e-001 Ad: 4.15e-001 Dobj: -1.7315108e+004 
Iter:  8 Ap: 3.54e-001 Pobj:  1.4419978e-001 Ad: 2.59e-001 Dobj: -3.0261606e+004 
Iter:  9 Ap: 2.06e-001 Pobj:  7.6159949e-002 Ad: 1.69e-001 Dobj: -5.1608031e+004 
Iter: 10 Ap: 1.76e-001 Pobj:  2.8225580e-002 Ad: 2.05e-001 Dobj: -1.0381123e+005 
Iter: 11 Ap: 4.38e-002 Pobj: -1.0884879e-001 Ad: 2.70e-002 Dobj: -1.5849860e+005 
Iter: 12 Ap: 7.41e-002 Pobj: -1.0017728e-001 Ad: 3.78e-002 Dobj: -2.2427085e+005 
Iter: 13 Ap: 8.04e-002 Pobj: -2.3143836e-001 Ad: 5.03e-002 Dobj: -3.1189763e+005 
Iter: 14 Ap: 3.12e-002 Pobj: -2.5591558e-001 Ad: 2.20e-002 Dobj: -3.9378658e+005 
Iter: 15 Ap: 5.86e-002 Pobj: -3.8882386e-001 Ad: 7.00e-002 Dobj: -7.6150734e+005 
Iter: 16 Ap: 2.04e-002 Pobj: -3.8027812e-001 Ad: 1.55e-002 Dobj: -8.7496163e+005 
Iter: 17 Ap: 6.65e-002 Pobj: -5.9772591e-001 Ad: 3.41e-002 Dobj: -1.1931748e+006 
Iter: 18 Ap: 3.03e-002 Pobj: -5.8196452e-001 Ad: 1.98e-002 Dobj: -1.4676507e+006 
Iter: 19 Ap: 5.93e-002 Pobj: -8.0068609e-001 Ad: 2.37e-002 Dobj: -1.8016859e+006 
Iter: 20 Ap: 7.13e-002 Pobj: -7.9588464e-001 Ad: 4.77e-002 Dobj: -2.3435909e+006 
Iter: 21 Ap: 5.84e-002 Pobj: -1.0256495e+000 Ad: 4.88e-002 Dobj: -2.9289293e+006 
Iter: 22 Ap: 6.33e-002 Pobj: -1.0293021e+000 Ad: 7.47e-002 Dobj: -4.2752813e+006 
Iter: 23 Ap: 2.00e-001 Pobj: -1.3404625e+000 Ad: 2.01e-001 Dobj: -4.5856521e+006 
Iter: 24 Ap: 2.42e-001 Pobj: -1.5124606e+000 Ad: 1.82e-001 Dobj: -4.5858794e+006 
Iter: 25 Ap: 3.10e-001 Pobj: -1.8383679e+000 Ad: 2.52e-001 Dobj: -4.4679781e+006 
Iter: 26 Ap: 3.52e-001 Pobj: -2.3632871e+000 Ad: 3.28e-001 Dobj: -4.1256414e+006 
Iter: 27 Ap: 2.59e-001 Pobj: -3.3773445e+000 Ad: 2.48e-001 Dobj: -3.8335879e+006 
Iter: 28 Ap: 2.04e-001 Pobj: -5.5014065e+000 Ad: 1.78e-001 Dobj: -3.5991974e+006 
Iter: 29 Ap: 1.31e-001 Pobj: -9.2674176e+000 Ad: 1.25e-001 Dobj: -3.4380223e+006 
Iter: 30 Ap: 8.62e-002 Pobj: -1.5659852e+001 Ad: 9.49e-002 Dobj: -3.3222368e+006 
Iter: 31 Ap: 5.82e-002 Pobj: -2.6072022e+001 Ad: 7.46e-002 Dobj: -3.2342994e+006 
Iter: 32 Ap: 4.05e-002 Pobj: -4.2447605e+001 Ad: 6.18e-002 Dobj: -3.1619966e+006 
Iter: 33 Ap: 2.97e-002 Pobj: -6.8060364e+001 Ad: 5.27e-002 Dobj: -3.0992810e+006 
Iter: 34 Ap: 2.29e-002 Pobj: -1.0870437e+002 Ad: 3.61e-002 Dobj: -3.0555436e+006 
Iter: 35 Ap: 1.94e-002 Pobj: -1.7519156e+002 Ad: 3.14e-002 Dobj: -3.0181803e+006 
Iter: 36 Ap: 2.71e-002 Pobj: -3.5585378e+002 Ad: 3.36e-002 Dobj: -2.9794972e+006 
Iter: 37 Ap: 3.16e-002 Pobj: -8.9136747e+002 Ad: 3.08e-002 Dobj: -2.9459964e+006 
Iter: 38 Ap: 2.98e-002 Pobj: -2.5515147e+003 Ad: 2.59e-002 Dobj: -2.9235276e+006 
Iter: 39 Ap: 3.32e-003 Pobj: -3.2749581e+003 Ad: 3.04e-002 Dobj: -2.9068907e+006 
Iter: 40 Ap: 1.71e-003 Pobj: -4.0293522e+003 Ad: 2.75e-002 Dobj: -2.8934632e+006 
Iter: 41 Ap: 7.13e-004 Pobj: -4.6852931e+003 Ad: 1.95e-002 Dobj: -2.8811600e+006 
Iter: 42 Ap: 6.28e-004 Pobj: -5.6563176e+003 Ad: 6.37e-003 Dobj: -2.8716478e+006 
Iter: 43 Ap: 7.59e-004 Pobj: -7.2852922e+003 Ad: 1.68e-002 Dobj: -2.8633124e+006 
Iter: 44 Ap: 5.97e-004 Pobj: -9.8923166e+003 Ad: 5.69e-003 Dobj: -2.8519349e+006 
Iter: 45 Ap: 2.75e-004 Pobj: -1.1609626e+004 Ad: 9.21e-003 Dobj: -2.8433925e+006 
Iter: 46 Ap: 7.91e-004 Pobj: -1.9123711e+004 Ad: 4.11e-003 Dobj: -2.8400379e+006 
Iter: 47 Ap: 1.49e-004 Pobj: -2.1788925e+004 Ad: 3.72e-003 Dobj: -2.8361820e+006 
Iter: 48 Ap: 1.91e-004 Pobj: -2.6028917e+004 Ad: 4.04e-003 Dobj: -2.8304991e+006 
Iter: 49 Ap: 4.92e-005 Pobj: -2.7432397e+004 Ad: 4.24e-003 Dobj: -2.8252920e+006 
Iter: 50 Ap: 3.18e-004 Pobj: -3.8306243e+004 Ad: 4.66e-003 Dobj: -2.8222317e+006 
Iter: 51 Ap: 4.26e-005 Pobj: -4.0857207e+004 Ad: 3.87e-003 Dobj: -2.8194880e+006 
Iter: 52 Ap: 2.99e-004 Pobj: -6.3688468e+004 Ad: 3.34e-003 Dobj: -2.8151265e+006 
Iter: 53 Ap: 1.67e-004 Pobj: -8.5640730e+004 Ad: 1.89e-003 Dobj: -2.8130403e+006 
Iter: 54 Ap: 1.51e-004 Pobj: -1.1470022e+005 Ad: 2.46e-003 Dobj: -2.8114466e+006 
Iter: 55 Ap: 1.61e-004 Pobj: -1.6358658e+005 Ad: 2.52e-003 Dobj: -2.8098889e+006 
Iter: 56 Ap: 1.05e-004 Pobj: -2.1860126e+005 Ad: 1.22e-003 Dobj: -2.8095144e+006 
Iter: 57 Ap: 3.61e-005 Pobj: -2.4668581e+005 Ad: 2.43e-003 Dobj: -2.8075831e+006 
Iter: 58 Ap: 7.26e-005 Pobj: -3.2183543e+005 Ad: 2.39e-003 Dobj: -2.8048472e+006 
Iter: 59 Ap: 3.28e-005 Pobj: -3.7297373e+005 Ad: 8.33e-004 Dobj: -2.8045873e+006 
Iter: 60 Ap: 9.91e-005 Pobj: -5.6430845e+005 Ad: 1.11e-003 Dobj: -2.8027903e+006 
Iter: 61 Ap: 5.65e-005 Pobj: -7.3405737e+005 Ad: 1.62e-003 Dobj: -2.8007675e+006 
Lack of progress.  Giving up!
Failure: return code is 7 
Primal objective value: -6.8060364e+001
Dual objective value: -3.0992810e+006
Relative primal infeasibility: 5.27e-002
Relative dual infeasibility: 3.71e+000
Real Relative Gap: -1.00e+000
XZ Relative Gap: 3.69e+000
DIMACS error measures: 7.57e-002 0.00e+000 6.40e+000 0.00e+000 -1.00e+000 3.69e+000
objective_value(model) = 68.06036356465847
68.06036356465847
```

## 2022.3.14

现在需要做的事情：
- 读CSDP的输出
- 检查约束生成的是不是正确（使用Mathematica）
- 拿nonlinear SDP的输出做benchmark

在`calculate-all-correlation-functions-in-xp-oscillator-for-benchmark-1.nb`中做benchmark。计算得到$\lang x^{2n} \rang$得到
```
0.301138, 0.253782, 0.343358, 0.598177, 1.31284, 3.40689, 9.75835, 
32.4658, 117.962, 451.727
```
$n$大的时候这个期望值大得离谱。无法理解。为谨慎起见最好同时做一下解方程解出来的$\lang x^{2n} \rang$。

以下是`stationary-schrodinger.jl`给出的列表
```
   0.3056965917452047
   0.2600508887447216
   0.3469360154545109
   0.6157259808259774
   1.3445134629688316
   3.4523744557836737
  10.122111626479073
  33.196322791312184
 119.94730115539195
 471.98535693014253
```
看起来是差不多的。

在`calculate-all-correlation-functions-in-xp-oscillator-for-benchmark-1.nb`中，运行
```
Table[comm[xntimes[i], H], {i, 1, 10} ]
```
能够得到
```
{2 I xpOpString[p], 2 xpOpString[] + 4 I xpOpString[x, p], 
 6 xpOpString[x] + 6 I xpOpString[x, x, p], 
 12 xpOpString[x, x] + 8 I xpOpString[x, x, x, p], 
 20 xpOpString[x, x, x] + 10 I xpOpString[x, x, x, x, p], 
 30 xpOpString[x, x, x, x] + 12 I xpOpString[x, x, x, x, x, p], 
 42 xpOpString[x, x, x, x, x] + 14 I xpOpString[x, x, x, x, x, x, p], 
 56 xpOpString[x, x, x, x, x, x] + 
  16 I xpOpString[x, x, x, x, x, x, x, p], 
 72 xpOpString[x, x, x, x, x, x, x] + 
  18 I xpOpString[x, x, x, x, x, x, x, x, p], 
 90 xpOpString[x, x, x, x, x, x, x, x] + 
  20 I xpOpString[x, x, x, x, x, x, x, x, x, p]}
```
然后可以看到，奇数个算符乘积的期望值是零这个事实是可以通过和哈密顿量的对易关系得到的。

感觉一种可能的做法是把所有的约束全部列出来，在`L_max`比较小的时候，这样能够判断什么约束是已经有了的。

一批含有一个$p$的算符的期望值：
```
{{xpOpString[x, p] -> 0. + 0.5 I, xpOpString[x, x, p] -> 0. + 0. I, 
  xpOpString[x, x, x, p] -> 0. + 0.451707 I, 
  xpOpString[x, x, x, x, p] -> 0. + 0. I, 
  xpOpString[x, x, x, x, x, p] -> 0. + 0.634454 I, 
  xpOpString[x, x, x, x, x, x, p] -> 0. + 0. I, 
  xpOpString[x, x, x, x, x, x, x, p] -> 0. + 1.20175 I, 
  xpOpString[x, x, x, x, x, x, x, x, p] -> 0. + 0. I, 
  xpOpString[x, x, x, x, x, x, x, x, x, p] -> 0. + 2.6918 I}}
```
等待和`stationary-schrodinger.jl`对比。

## 2022.3.15

对了这里头还有一个点，就是一维束缚态波函数总是实的，而$p = - \ii \grad$，所以含有奇数个$p$的算符序列的期望值是虚的，含有偶数个$p$的算符序列的期望值是实的。

找文献找到这个：http://www.optimization-online.org/DB_FILE/2002/10/551.pdf

lacking of progress是独立于max iteration times的。不知道是不是和初始值有关系？？？

想起来还有一个事情。~被我这么一弄好像`M`没法保证是一个对称的矩阵啊~这个没有问题因为`j`从`i`开始。在一次失败的优化之后运行如下诊断代码：
```julia
Mcons = Matrix(undef, 2 * (L_max + 1)^2, 2 * (L_max + 1)^2)
for i in 1 : (L_max + 1)^2
    for j in i : (L_max + 1)^2
        op1_idx = M_index_to_xpopstr_index[i]
        op2_idx = M_index_to_xpopstr_index[j]
        op1_idx_xpower = index_to_xpower(op1_idx)
        op1_idx_ppower = index_to_ppower(op1_idx)
        op2_idx_xpower = index_to_xpower(op2_idx)
        op2_idx_ppower = index_to_ppower(op2_idx)
        op_ij = xpopstr_normal_ord(op1_idx_xpower, op1_idx_ppower, op2_idx_xpower, op2_idx_ppower)

        real_part = transpose(real(op_ij)) * xpopstr_basis_real * I22
        imag_part = transpose(imag(op_ij)) * xpopstr_basis_imag * Im22
        Mcons[2i - 1 : 2i, 2j - 1 : 2j] = real_part + imag_part
    end
end
Mcons
```
马上发现两个大问题，~一个是有一些矩阵元好像没有被遍历到~（这个是正确的，这样直接解决了前述不对称的问题），一个是
```julia
Mcons[1:2, 3:4]
```
给出
```
2×2 Array{Any,2}:
 xpopstr_expected[1]  0
 0                    xpopstr_expected[1]
```
虚部去哪儿了？？？

我在这里大概率还是犯了虚部实部分开来算的问题，交叉项被忽略了。如下诊断代码
```julia
Mcons = Matrix(undef, 2 * (L_max + 1)^2, 2 * (L_max + 1)^2)
for i in 1 : (L_max + 1)^2
    for j in i : (L_max + 1)^2
        op1_idx = M_index_to_xpopstr_index[i]
        op2_idx = M_index_to_xpopstr_index[j]
        op1_idx_xpower = index_to_xpower(op1_idx)
        op1_idx_ppower = index_to_ppower(op1_idx)
        op2_idx_xpower = index_to_xpower(op2_idx)
        op2_idx_ppower = index_to_ppower(op2_idx)
        op_ij = xpopstr_normal_ord(op1_idx_xpower, op1_idx_ppower, op2_idx_xpower, op2_idx_ppower)

        Mcons[2i - 1 : 2i, 2j - 1 : 2j] = complex_to_mat(op_ij)
    end
end
Mcons
```
至少能够给出还算正常的结果。将其复制回`jump-oscillator-2.jl`。

收敛是收敛了，但是结果明显不正确，估计还是有bug。下面要怎么检查呢……首先尝试把`L_max`加大，至少看看有没有特别“硬”的错误。
这是`L_max = 8`的结果：
```
CSDP 6.2.0
Iter:  0 Ap: 0.00e+000 Pobj:  0.0000000e+000 Ad: 0.00e+000 Dobj:  0.0000000e+000 
Iter:  1 Ap: 3.53e-002 Pobj: -1.1371252e-001 Ad: 3.32e-002 Dobj: -1.2633893e+010 
Iter:  2 Ap: 3.86e-002 Pobj:  1.9771642e-001 Ad: 4.07e-002 Dobj: -3.1670609e+010 
Iter:  3 Ap: 9.01e-003 Pobj:  5.9339128e-001 Ad: 1.08e-002 Dobj: -1.0581058e+011 
Iter:  4 Ap: 8.15e-003 Pobj:  6.9908274e-001 Ad: 3.86e-003 Dobj: -1.3540563e+011 
Iter:  5 Ap: 8.69e-003 Pobj:  1.0760270e+000 Ad: 5.33e-003 Dobj: -1.5910648e+011 
Iter:  6 Ap: 2.85e-002 Pobj:  7.1915772e-001 Ad: 1.61e-002 Dobj: -1.7208901e+011 
Iter:  7 Ap: 3.75e-003 Pobj:  7.1225932e-001 Ad: 2.82e-003 Dobj: -1.9343672e+011 
Iter:  8 Ap: 5.12e-003 Pobj:  5.8690978e-001 Ad: 6.62e-003 Dobj: -2.4414726e+011 
Iter:  9 Ap: 6.00e-003 Pobj:  4.6037927e-001 Ad: 7.88e-003 Dobj: -2.7921638e+011 
Iter: 10 Ap: 4.93e-003 Pobj:  3.3893493e-001 Ad: 9.15e-003 Dobj: -3.3357332e+011 
Iter: 11 Ap: 5.72e-002 Pobj:  2.5671945e-001 Ad: 4.12e-002 Dobj: -4.8430713e+011 
Iter: 12 Ap: 5.19e-003 Pobj:  3.0155970e-001 Ad: 6.47e-003 Dobj: -5.3597790e+011 
Iter: 13 Ap: 6.73e-003 Pobj:  5.3901428e-001 Ad: 2.88e-003 Dobj: -5.4167460e+011 
Iter: 14 Ap: 9.73e-003 Pobj:  4.8494183e-001 Ad: 7.08e-003 Dobj: -5.7151539e+011 
Iter: 15 Ap: 2.24e-002 Pobj:  9.0724511e-001 Ad: 1.61e-002 Dobj: -6.6037626e+011 
Iter: 16 Ap: 1.21e-002 Pobj:  9.1920960e-001 Ad: 2.53e-002 Dobj: -7.8849577e+011 
Iter: 17 Ap: 2.93e-002 Pobj:  1.1340522e+000 Ad: 9.34e-003 Dobj: -8.3302010e+011 
Iter: 18 Ap: 1.46e-002 Pobj:  1.1260562e+000 Ad: 1.97e-002 Dobj: -9.1213945e+011 
Iter: 19 Ap: 6.42e-002 Pobj:  9.7611595e-001 Ad: 6.37e-002 Dobj: -1.1628454e+012 
Iter: 20 Ap: 5.35e-003 Pobj:  8.2967511e-001 Ad: 1.26e-002 Dobj: -1.2033384e+012 
Iter: 21 Ap: 4.83e-002 Pobj:  6.5228258e-001 Ad: 1.86e-002 Dobj: -1.2662103e+012 
Iter: 22 Ap: 2.72e-002 Pobj:  6.0957246e-001 Ad: 8.89e-003 Dobj: -1.2920288e+012 
Iter: 23 Ap: 7.35e-002 Pobj:  7.4214782e-001 Ad: 4.21e-002 Dobj: -1.4110771e+012 
Iter: 24 Ap: 7.01e-003 Pobj:  7.5657707e-001 Ad: 8.20e-003 Dobj: -1.4182229e+012 
Iter: 25 Ap: 5.33e-002 Pobj:  1.3258253e+000 Ad: 1.07e-002 Dobj: -1.4466905e+012 
Iter: 26 Ap: 3.02e-002 Pobj:  1.0936977e+000 Ad: 5.41e-002 Dobj: -1.5520664e+012 
Iter: 27 Ap: 3.86e-002 Pobj:  9.7203783e-001 Ad: 4.10e-002 Dobj: -1.6219063e+012 
Iter: 28 Ap: 9.25e-002 Pobj:  1.2646083e+000 Ad: 4.03e-002 Dobj: -1.6847849e+012 
Iter: 29 Ap: 7.59e-002 Pobj:  8.8542941e-001 Ad: 6.79e-002 Dobj: -1.7722829e+012 
Iter: 30 Ap: 1.19e-001 Pobj:  1.1485405e+000 Ad: 1.03e-001 Dobj: -1.8809732e+012 
Iter: 31 Ap: 2.22e-002 Pobj:  1.1202238e+000 Ad: 3.66e-002 Dobj: -1.9048941e+012 
Iter: 32 Ap: 1.00e-001 Pobj:  1.0405420e+000 Ad: 1.25e-001 Dobj: -1.9847845e+012 
Iter: 33 Ap: 5.72e-002 Pobj:  9.8523370e-001 Ad: 6.35e-002 Dobj: -2.0160262e+012 
Iter: 34 Ap: 1.64e-001 Pobj:  8.7117162e-001 Ad: 9.88e-002 Dobj: -2.0593502e+012 
Iter: 35 Ap: 1.60e-001 Pobj:  9.3937496e-001 Ad: 1.51e-001 Dobj: -2.1116801e+012 
Iter: 36 Ap: 8.42e-002 Pobj:  8.5832673e-001 Ad: 1.20e-001 Dobj: -2.1404004e+012 
Iter: 37 Ap: 2.40e-001 Pobj:  8.8393094e-001 Ad: 1.05e-001 Dobj: -2.1563065e+012 
Iter: 38 Ap: 1.18e-001 Pobj:  7.5037751e-001 Ad: 9.42e-002 Dobj: -2.1622270e+012 
Iter: 39 Ap: 1.89e-001 Pobj:  8.5428144e-001 Ad: 2.02e-001 Dobj: -2.1653518e+012 
Iter: 40 Ap: 3.69e-001 Pobj:  5.8646190e-001 Ad: 1.51e-001 Dobj: -2.1703108e+012 
Iter: 41 Ap: 1.39e-001 Pobj:  6.8988975e-001 Ad: 1.92e-001 Dobj: -2.1746625e+012 
Iter: 42 Ap: 1.24e-001 Pobj:  8.2613044e-001 Ad: 1.09e-001 Dobj: -2.1759579e+012 
Iter: 43 Ap: 2.69e-001 Pobj:  7.8100792e-001 Ad: 1.44e-001 Dobj: -2.1783869e+012 
Iter: 44 Ap: 2.42e-001 Pobj:  6.3088385e-001 Ad: 2.78e-001 Dobj: -2.1839801e+012 
Iter: 45 Ap: 2.24e-001 Pobj:  5.5384965e-001 Ad: 2.24e-001 Dobj: -2.1883305e+012 
Iter: 46 Ap: 1.93e-001 Pobj:  4.2565861e-001 Ad: 2.26e-001 Dobj: -2.1933863e+012 
Iter: 47 Ap: 1.73e-001 Pobj:  3.3185990e-001 Ad: 1.39e-001 Dobj: -2.1972666e+012 
Iter: 48 Ap: 3.42e-001 Pobj:  1.7940275e-001 Ad: 2.89e-001 Dobj: -2.2055049e+012 
Iter: 49 Ap: 2.40e-001 Pobj:  2.0283503e-001 Ad: 2.36e-001 Dobj: -2.2037523e+012 
Iter: 50 Ap: 2.77e-001 Pobj:  1.8111233e-001 Ad: 1.91e-001 Dobj: -2.2009851e+012 
Iter: 51 Ap: 2.97e-001 Pobj:  1.1279881e-001 Ad: 2.24e-001 Dobj: -2.1962655e+012 
Iter: 52 Ap: 2.96e-001 Pobj:  1.1546111e-001 Ad: 2.79e-001 Dobj: -2.1916748e+012 
Iter: 53 Ap: 2.02e-001 Pobj:  1.2508751e-001 Ad: 1.97e-001 Dobj: -2.1875736e+012 
Iter: 54 Ap: 3.24e-001 Pobj:  1.3121922e-001 Ad: 3.22e-001 Dobj: -2.1845570e+012 
Iter: 55 Ap: 2.36e-001 Pobj:  1.1395400e-001 Ad: 3.02e-001 Dobj: -2.1827091e+012 
Iter: 56 Ap: 2.85e-001 Pobj:  8.6176057e-002 Ad: 2.76e-001 Dobj: -2.1837005e+012 
Iter: 57 Ap: 2.91e-001 Pobj:  8.6935895e-002 Ad: 2.46e-001 Dobj: -2.1859433e+012 
Iter: 58 Ap: 3.76e-001 Pobj:  1.6136500e-002 Ad: 3.31e-001 Dobj: -2.1879369e+012 
Iter: 59 Ap: 4.08e-001 Pobj:  4.7631888e-002 Ad: 4.31e-001 Dobj: -2.1896275e+012 
Declaring primal infeasibility.
Success: SDP is primal infeasible
Certificate of primal infeasibility: a'*y=-1.00000e+000, ||A'(y)-Z||=8.04969e-009
objective_value(model) = -0.04763188796641771
```
因此是可以收敛的，没有诸如指标越界之类的问题。

下一步要做的事情是观察到底何处出了差错。检查各个关联函数的值。运行诊断代码
```julia
for i in 1 : L_max
    println("i = $i : $(value(xpopstr_expected_real_imag_parts(xpopstr_index(i, 0), :real)))")
end
```
得到
```
i = 1 : 0.0
i = 2 : -0.022826707994681783
i = 3 : 1.4551915228366852e-11
i = 4 : -0.0067335206404095516
i = 5 : -1.4551915228366852e-11
i = 6 : -0.17618527909507975
i = 7 : 0.0
i = 8 : -0.33923411423165817
```
这里`i = `给出$\lang x^i \rang$。没有什么要求$\lang x^2 \rang$大于零吗？

要么是我程序写错了，什么约束不对，要么是程序实际上是对的，但是`L_max`不够大。

将`L_max`设置成`10`差点让电脑炸掉。还是拿到服务器上去算吧……

创建`jump-oscillator-2-cluster-version.jl`，它是`jump-oscillator-2.jl`的复制品。

~by the way:好像CSDP的优化还不是确定性的。运行一遍`jump-oscillator-2-cluster-version.jl`，给出~
草，怎么会"SDP is primal infeasible"

以下是输出：
```
CSDP 6.2.0
Iter:  0 Ap: 0.00e+000 Pobj:  0.0000000e+000 Ad: 0.00e+000 Dobj:  0.0000000e+000
Iter:  1 Ap: 3.53e-002 Pobj: -1.1371252e-001 Ad: 3.32e-002 Dobj: -1.2633893e+010
Iter:  2 Ap: 3.86e-002 Pobj:  1.9771642e-001 Ad: 4.07e-002 Dobj: -3.1670609e+010
Iter:  3 Ap: 9.01e-003 Pobj:  5.9339128e-001 Ad: 1.08e-002 Dobj: -1.0581058e+011
Iter:  4 Ap: 8.15e-003 Pobj:  6.9908274e-001 Ad: 3.86e-003 Dobj: -1.3540563e+011
Iter:  5 Ap: 8.69e-003 Pobj:  1.0760270e+000 Ad: 5.33e-003 Dobj: -1.5910648e+011
Iter:  6 Ap: 2.85e-002 Pobj:  7.1915772e-001 Ad: 1.61e-002 Dobj: -1.7208901e+011
Iter:  7 Ap: 3.75e-003 Pobj:  7.1225932e-001 Ad: 2.82e-003 Dobj: -1.9343672e+011 
Iter:  8 Ap: 5.12e-003 Pobj:  5.8690978e-001 Ad: 6.62e-003 Dobj: -2.4414726e+011 
Iter:  9 Ap: 6.04e-003 Pobj:  4.5492560e-001 Ad: 8.26e-003 Dobj: -2.7346134e+011 
Iter: 10 Ap: 5.53e-003 Pobj:  3.4617285e-001 Ad: 1.06e-002 Dobj: -3.4826094e+011 
Iter: 11 Ap: 1.87e-002 Pobj:  2.5723368e-001 Ad: 1.77e-002 Dobj: -4.4285256e+011 
Iter: 12 Ap: 1.69e-002 Pobj:  3.7392408e-001 Ad: 1.14e-002 Dobj: -5.3041260e+011 
Iter: 13 Ap: 5.78e-003 Pobj:  4.8892867e-001 Ad: 2.18e-003 Dobj: -5.3740973e+011 
Iter: 14 Ap: 2.58e-002 Pobj:  5.4932685e-001 Ad: 1.58e-002 Dobj: -6.0829315e+011 
Iter: 15 Ap: 1.25e-002 Pobj:  7.4332223e-001 Ad: 9.21e-003 Dobj: -6.5925436e+011 
Iter: 16 Ap: 2.58e-002 Pobj:  1.0659186e+000 Ad: 1.47e-002 Dobj: -7.2977706e+011 
Iter: 17 Ap: 7.13e-003 Pobj:  9.3901698e-001 Ad: 1.77e-002 Dobj: -8.1382734e+011 
Iter: 18 Ap: 2.46e-002 Pobj:  1.0895576e+000 Ad: 1.86e-002 Dobj: -8.9799562e+011 
Iter: 19 Ap: 3.17e-002 Pobj:  1.1073767e+000 Ad: 2.27e-002 Dobj: -9.8638733e+011 
Iter: 20 Ap: 6.01e-002 Pobj:  1.0032226e+000 Ad: 4.85e-002 Dobj: -1.1719769e+012 
Iter: 21 Ap: 1.01e-002 Pobj:  8.6941308e-001 Ad: 1.04e-002 Dobj: -1.2008798e+012 
Iter: 22 Ap: 4.73e-002 Pobj:  8.2370245e-001 Ad: 6.05e-002 Dobj: -1.3887016e+012 
Iter: 23 Ap: 6.23e-002 Pobj:  8.2167299e-001 Ad: 2.39e-002 Dobj: -1.4551220e+012 
Iter: 24 Ap: 9.18e-003 Pobj:  8.6537322e-001 Ad: 1.17e-002 Dobj: -1.4673988e+012 
Iter: 25 Ap: 5.38e-002 Pobj:  1.1133512e+000 Ad: 7.67e-003 Dobj: -1.4866725e+012 
Iter: 26 Ap: 3.33e-002 Pobj:  1.1622298e+000 Ad: 1.67e-002 Dobj: -1.5208576e+012 
Iter: 27 Ap: 6.70e-002 Pobj:  1.1839656e+000 Ad: 7.12e-002 Dobj: -1.6466306e+012 
Iter: 28 Ap: 9.10e-002 Pobj:  1.0550776e+000 Ad: 7.59e-002 Dobj: -1.7556176e+012 
Iter: 29 Ap: 8.95e-002 Pobj:  1.0942434e+000 Ad: 8.31e-002 Dobj: -1.8509235e+012 
Iter: 30 Ap: 4.05e-003 Pobj:  1.0935390e+000 Ad: 5.87e-002 Dobj: -1.8969044e+012 
Iter: 31 Ap: 1.70e-001 Pobj:  1.0790382e+000 Ad: 5.64e-002 Dobj: -1.9443992e+012 
Iter: 32 Ap: 4.11e-002 Pobj:  1.0561925e+000 Ad: 6.42e-002 Dobj: -1.9744064e+012 
Iter: 33 Ap: 1.51e-001 Pobj:  9.6686382e-001 Ad: 8.25e-002 Dobj: -2.0166466e+012 
Iter: 34 Ap: 1.04e-001 Pobj:  8.9213675e-001 Ad: 1.06e-001 Dobj: -2.0660737e+012 
Iter: 35 Ap: 1.43e-001 Pobj:  9.3868120e-001 Ad: 1.06e-001 Dobj: -2.0971558e+012 
Iter: 36 Ap: 1.55e-001 Pobj:  7.5889202e-001 Ad: 1.51e-001 Dobj: -2.1358388e+012 
Iter: 37 Ap: 2.70e-001 Pobj:  6.8179533e-001 Ad: 9.72e-002 Dobj: -2.1524306e+012 
Iter: 38 Ap: 1.53e-001 Pobj:  6.7738180e-001 Ad: 1.12e-001 Dobj: -2.1556973e+012 
Iter: 39 Ap: 1.90e-001 Pobj:  7.6388764e-001 Ad: 2.53e-001 Dobj: -2.1628003e+012 
Iter: 40 Ap: 3.50e-001 Pobj:  6.1602193e-001 Ad: 1.46e-001 Dobj: -2.1678800e+012 
Iter: 41 Ap: 2.46e-001 Pobj:  6.1710460e-001 Ad: 2.31e-001 Dobj: -2.1743039e+012 
Iter: 42 Ap: 1.86e-001 Pobj:  6.8689723e-001 Ad: 1.40e-001 Dobj: -2.1764969e+012 
Iter: 43 Ap: 3.14e-001 Pobj:  4.5795904e-001 Ad: 1.70e-001 Dobj: -2.1803305e+012 
Iter: 44 Ap: 1.93e-001 Pobj:  3.8764601e-001 Ad: 2.46e-001 Dobj: -2.1876301e+012 
Iter: 45 Ap: 1.86e-001 Pobj:  3.3430853e-001 Ad: 1.62e-001 Dobj: -2.1937029e+012 
Iter: 46 Ap: 2.74e-001 Pobj:  1.7238648e-001 Ad: 2.57e-001 Dobj: -2.2037269e+012 
Iter: 47 Ap: 2.62e-001 Pobj:  1.6717003e-001 Ad: 1.83e-001 Dobj: -2.2100008e+012 
Iter: 48 Ap: 3.22e-001 Pobj:  1.5063310e-001 Ad: 1.95e-001 Dobj: -2.2101090e+012 
Iter: 49 Ap: 2.89e-001 Pobj:  1.0433502e-001 Ad: 2.05e-001 Dobj: -2.2071497e+012
Iter: 50 Ap: 2.83e-001 Pobj:  7.3228531e-002 Ad: 1.73e-001 Dobj: -2.2021428e+012
Iter: 51 Ap: 3.02e-001 Pobj:  8.2581573e-002 Ad: 2.75e-001 Dobj: -2.1968855e+012
Iter: 52 Ap: 2.12e-001 Pobj:  4.7739664e-002 Ad: 1.44e-001 Dobj: -2.1931547e+012
Iter: 53 Ap: 2.63e-001 Pobj:  3.2399536e-002 Ad: 2.63e-001 Dobj: -2.1894749e+012
Iter: 54 Ap: 2.59e-001 Pobj:  5.1823445e-002 Ad: 3.38e-001 Dobj: -2.1869945e+012
Iter: 55 Ap: 2.82e-001 Pobj:  1.7404358e-002 Ad: 3.18e-001 Dobj: -2.1910413e+012
Iter: 56 Ap: 2.83e-001 Pobj:  3.7743292e-002 Ad: 2.73e-001 Dobj: -2.2040088e+012
Iter: 57 Ap: 4.26e-001 Pobj:  2.6458428e-002 Ad: 3.61e-001 Dobj: -2.2079154e+012
Iter: 58 Ap: 4.25e-001 Pobj:  2.8240029e-002 Ad: 4.34e-001 Dobj: -2.2110354e+012
Iter: 59 Ap: 4.35e-001 Pobj:  8.1989347e-002 Ad: 3.77e-001 Dobj: -2.2112856e+012
Declaring primal infeasibility.
Success: SDP is primal infeasible
Certificate of primal infeasibility: a'*y=-1.00000e+000, ||A'(y)-Z||=6.46771e-009
objective_value(model) = -0.08198934656684287
```

这么说估计是哪里写错了。

## 2022.3.16

即使将`L_max`设置为5，同样有`Declaring primal infeasibility.`。

先看看这个时候约束对不对吧。运行诊断代码
```julia
for i in 1 : L_max
    println("i = $i : $(value(xpopstr_expected_real_imag_parts(xpopstr_index(i, 0), :real))) + $(value(xpopstr_expected_real_imag_parts(xpopstr_index(i, 0), :imag))) im ")
end
```
得到
```
i = 1 : -1.7763568394002505e-15 + 0.0 im 
i = 2 : -0.05054150351097597 + -3.3311966518567715e-9 im
i = 3 : 0.0 + 0.0 im
i = 4 : 0.0558192989219819 + 1.2002825400259098e-9 im
i = 5 : 0.0 + 0.0 im
```
因此奇数$\lang x^n \rang$为零的约束施加得很好，且$\Im \lang x^n \rang$几乎为零。

我们再来看$\land x^n p \rang$。运行诊断代码
```julia
for i in 1 : L_max
    println("i = $i : $(value(xpopstr_expected_real_imag_parts(xpopstr_index(i, 1), :real))) + $(value(xpopstr_expected_real_imag_parts(xpopstr_index(i, 1), :imag))) im ")
end
```
得到
```
i = 1 : 0.0743289004671901 + -0.38974764309131515 im 
i = 2 : 0.0 + 0.0 im
i = 3 : -0.07827649436848816 + 0.07794952716507275 im
i = 4 : 0.0 + 0.0 im
i = 5 : 0.1423146029647775 + 0.008878796747188389 im
```
说明约束$\Re \lang x^{2n-1}p \rang = 0$并没有得到实施。另一方面，$\lang x^{2n} p \rang = 0$得到了实施。
前者是基态才有的性质，后者是哈密顿量决定的，因此至少一部分由哈密顿量施加的约束得到了实施，但是优化存在难度。

为了分析$M$矩阵是否正确构造，执行诊断代码
```julia
storage_name = "D:\\Projects\\numerical-boostrap\\oscillator-simple-prototype\\"
open(storage_name * "M-mat-element-oscillator-2-cluster-version-2022-3-16.txt", "w") do file
    for i in 1 : (L_max + 1)^2
        for j in i : (L_max + 1)^2
            op1_idx = M_index_to_xpopstr_index[i]
            op2_idx = M_index_to_xpopstr_index[j]
            op1_idx_xpower = index_to_xpower(op1_idx)
            op1_idx_ppower = index_to_ppower(op1_idx)
            op2_idx_xpower = index_to_xpower(op2_idx)
            op2_idx_ppower = index_to_ppower(op2_idx)
            op_ij = xpopstr_normal_ord(op1_idx_xpower, op1_idx_ppower, op2_idx_xpower, op2_idx_ppower)
    
            println(file, 
                "i = $i  j = $j  x^$op1_idx_xpower p^$op1_idx_ppower x^$op2_idx_xpower p^$op2_idx_ppower ")
            println(file, complex_to_mat(op_ij)[1, 1])
            println(file, complex_to_mat(op_ij)[1, 2])
            println(file, complex_to_mat(op_ij)[2, 1])
            println(file, complex_to_mat(op_ij)[2, 2])
        end
    end
end
```
所得结果保存于`M-mat-element-oscillator-2-cluster-version-2022-3-16.txt`中。

为了分析是否$M$矩阵被正确的建立了，如下做抽检：从`M-mat-element-oscillator-2-cluster-version-2022-3-16.txt`中取样
```
i = 24  j = 36  x^3 p^5 x^5 p^5 
-120 xpopstr_expected[75] + 600 xpopstr_expected[99] + 600 xpopstr_expected[123] - 200 xpopstr_expected[147] - 25 xpopstr_expected[171] + xpopstr_expected[195]
120 xpopstr_expected[76] - 600 xpopstr_expected[100] - 600 xpopstr_expected[124] + 200 xpopstr_expected[148] + 25 xpopstr_expected[172] - xpopstr_expected[196]
-120 xpopstr_expected[76] + 600 xpopstr_expected[100] + 600 xpopstr_expected[124] - 200 xpopstr_expected[148] - 25 xpopstr_expected[172] + xpopstr_expected[196]
-120 xpopstr_expected[75] + 600 xpopstr_expected[99] + 600 xpopstr_expected[123] - 200 xpopstr_expected[147] - 25 xpopstr_expected[171] + xpopstr_expected[195]
```
这里涉及到的`xpopstr_expected`变量的index（从1开始）包括

75, 99, 123, 147, 171, 195

使用如下诊断代码将它们转化成$x$和$p$的幂次：
```julia
xpopstr_expected_indices = [75, 99, 123, 147, 171, 195]
for idx in xpopstr_expected_indices
    x_power = index_to_xpower(Int((idx + 1) / 2))
    p_power = index_to_ppower(Int((idx + 1) / 2))
    println("$idx  x^$x_power p^$p_power")
end
```
得到
```
75  x^3 p^5
99  x^4 p^6
123  x^5 p^7
147  x^6 p^8
171  x^7 p^9
195  x^8 p^10
```
解析计算（见`2022-3-16.nb`）给出
$$
\left(
\begin{array}{cc}
 120 \Im\left(x^3.p^5\right)-600 \Im\left(x^5.p^7\right)+25 \Im\left(x^7.p^9\right)+600
   \Re\left(x^4.p^6\right)-200 \Re\left(x^6.p^8\right)+\Re\left(x^8.p^{10}\right) & -600
   \Im\left(x^4.p^6\right)+200 \Im\left(x^6.p^8\right)-\Im\left(x^8.p^{10}\right)+120
   \Re\left(x^3.p^5\right)-600 \Re\left(x^5.p^7\right)+25 \Re\left(x^7.p^9\right) \\
 600 \Im\left(x^4.p^6\right)-200 \Im\left(x^6.p^8\right)+\Im\left(x^8.p^{10}\right)-120
   \Re\left(x^3.p^5\right)+600 \Re\left(x^5.p^7\right)-25 \Re\left(x^7.p^9\right) & 120
   \Im\left(x^3.p^5\right)-600 \Im\left(x^5.p^7\right)+25 \Im\left(x^7.p^9\right)+600
   \Re\left(x^4.p^6\right)-200 \Re\left(x^6.p^8\right)+\Re\left(x^8.p^{10}\right) \\
\end{array}
\right)
$$
卧槽好像还真的对不上……