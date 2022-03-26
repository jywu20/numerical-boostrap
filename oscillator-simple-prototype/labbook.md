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

结果发现如下代码：
```julia
function complex_to_mat(coefficients)
    real_part = transpose(real(coefficients))
    imag_part = transpose(imag(coefficients))
    real_part_mat_version = map(x -> x * I22, real_part)
    imag_part_mat_version = map(x -> x * I22, imag_part)
    (real_part_mat_version + imag_part_mat_version) * (xpopstr_basis_real + xpopstr_basis_imag)
end
```
中第二行的`I22`不对，应该是`Im22`。

然而修改之后仍然提示`Declaring primal infeasibility.`……

重复诊断
```julia
storage_name = "D:\\Projects\\numerical-boostrap\\oscillator-simple-prototype\\"
open(storage_name * "M-mat-element-oscillator-2-cluster-version-2022-3-16-2.txt", "w") do file
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
取样：
```
i = 24  j = 36  x^3 p^5 x^5 p^5 
120 xpopstr_expected[76] + 600 xpopstr_expected[99] - 600 xpopstr_expected[124] - 200 xpopstr_expected[147] + 25 xpopstr_expected[172] + xpopstr_expected[195]
120 xpopstr_expected[75] - 600 xpopstr_expected[100] - 600 xpopstr_expected[123] + 200 xpopstr_expected[148] + 25 xpopstr_expected[171] - xpopstr_expected[196]
-120 xpopstr_expected[75] + 600 xpopstr_expected[100] + 600 xpopstr_expected[123] - 200 xpopstr_expected[148] - 25 xpopstr_expected[171] + xpopstr_expected[196]
120 xpopstr_expected[76] + 600 xpopstr_expected[99] - 600 xpopstr_expected[124] - 200 xpopstr_expected[147] + 25 xpopstr_expected[172] + xpopstr_expected[195]
```
这回似乎是正确的。

整体运行`jump-oscillator-2-cluster-version.jl`，仍然存在“infeasible”的输出：
```
CSDP 6.2.0
Iter:  0 Ap: 0.00e+000 Pobj:  0.0000000e+000 Ad: 0.00e+000 Dobj:  0.0000000e+000 
Iter:  1 Ap: 8.18e-002 Pobj: -1.3962147e-002 Ad: 9.43e-002 Dobj: -2.1809110e+010 
Iter:  2 Ap: 1.36e-001 Pobj: -4.2888011e-001 Ad: 6.88e-002 Dobj: -3.3221426e+010 
Iter:  3 Ap: 2.81e-002 Pobj: -1.0940146e+000 Ad: 2.17e-002 Dobj: -7.5508937e+010 
Iter:  4 Ap: 3.81e-002 Pobj: -8.7152731e-001 Ad: 7.41e-002 Dobj: -2.3816175e+011 
Iter:  5 Ap: 6.80e-002 Pobj:  1.4322870e+000 Ad: 5.05e-002 Dobj: -3.3485435e+011 
Iter:  6 Ap: 5.87e-002 Pobj:  3.5434138e+000 Ad: 6.03e-002 Dobj: -4.0310930e+011 
Iter:  7 Ap: 5.16e-002 Pobj:  3.6233577e+000 Ad: 4.89e-002 Dobj: -4.4682471e+011 
Iter:  8 Ap: 1.08e-001 Pobj:  4.1490790e+000 Ad: 9.53e-002 Dobj: -5.1735824e+011 
Iter:  9 Ap: 9.35e-002 Pobj:  2.4398432e+000 Ad: 8.56e-002 Dobj: -5.6546517e+011 
Iter: 10 Ap: 4.92e-002 Pobj:  2.4157495e+000 Ad: 5.79e-002 Dobj: -5.8655939e+011 
Iter: 11 Ap: 2.17e-001 Pobj:  1.3680341e+000 Ad: 1.11e-001 Dobj: -6.2108891e+011 
Iter: 12 Ap: 1.06e-001 Pobj:  1.6787347e+000 Ad: 1.37e-001 Dobj: -6.3750707e+011 
Iter: 13 Ap: 1.95e-001 Pobj:  1.5151120e+000 Ad: 9.04e-002 Dobj: -6.4720303e+011 
Iter: 14 Ap: 1.52e-001 Pobj:  3.3565862e-001 Ad: 1.63e-001 Dobj: -6.6165588e+011 
Iter: 15 Ap: 1.58e-001 Pobj:  7.6631334e-002 Ad: 1.44e-001 Dobj: -6.7205402e+011 
Iter: 16 Ap: 2.63e-001 Pobj: -3.0307305e-001 Ad: 1.90e-001 Dobj: -6.8199187e+011 
Iter: 17 Ap: 2.26e-001 Pobj: -8.9357966e-001 Ad: 1.75e-001 Dobj: -6.8839332e+011 
Iter: 18 Ap: 2.10e-001 Pobj: -6.7280867e-001 Ad: 1.11e-001 Dobj: -6.9183384e+011 
Iter: 19 Ap: 1.92e-001 Pobj: -8.8180127e-001 Ad: 1.65e-001 Dobj: -6.9574259e+011 
Iter: 20 Ap: 2.19e-001 Pobj: -8.0838883e-001 Ad: 7.28e-002 Dobj: -6.9709513e+011 
Iter: 21 Ap: 2.02e-001 Pobj: -6.8410885e-001 Ad: 2.39e-001 Dobj: -7.0181048e+011 
Iter: 22 Ap: 1.82e-001 Pobj: -7.1210439e-001 Ad: 1.69e-001 Dobj: -7.0451834e+011 
Iter: 23 Ap: 3.75e-001 Pobj: -7.1955460e-001 Ad: 2.81e-001 Dobj: -7.0748138e+011 
Iter: 24 Ap: 2.36e-001 Pobj: -3.4911901e-001 Ad: 1.95e-001 Dobj: -7.0872672e+011 
Iter: 25 Ap: 2.04e-001 Pobj: -2.5080693e-001 Ad: 1.98e-001 Dobj: -7.0956273e+011 
Iter: 26 Ap: 2.98e-001 Pobj: -1.5662899e-001 Ad: 1.83e-001 Dobj: -7.1041007e+011 
Iter: 27 Ap: 2.47e-001 Pobj:  7.5027960e-004 Ad: 1.82e-001 Dobj: -7.1092184e+011 
Iter: 28 Ap: 2.80e-001 Pobj:  2.6233834e-002 Ad: 1.97e-001 Dobj: -7.1127606e+011 
Iter: 29 Ap: 3.27e-001 Pobj:  1.4558903e-003 Ad: 2.53e-001 Dobj: -7.1143068e+011 
Iter: 30 Ap: 3.14e-001 Pobj: -1.8949497e-002 Ad: 1.87e-001 Dobj: -7.1149932e+011 
Iter: 31 Ap: 3.63e-001 Pobj: -1.4733383e-002 Ad: 2.84e-001 Dobj: -7.1157593e+011 
Iter: 32 Ap: 2.69e-001 Pobj:  2.1126149e-003 Ad: 2.37e-001 Dobj: -7.1164616e+011 
Iter: 33 Ap: 2.84e-001 Pobj: -1.0628196e-001 Ad: 2.08e-001 Dobj: -7.1179961e+011 
Iter: 34 Ap: 2.24e-001 Pobj: -1.8582056e-001 Ad: 1.78e-001 Dobj: -7.1184138e+011 
Iter: 35 Ap: 1.22e-001 Pobj: -1.8141766e-001 Ad: 1.07e-001 Dobj: -7.1186888e+011 
Iter: 36 Ap: 1.90e-001 Pobj: -1.6999011e-001 Ad: 8.87e-002 Dobj: -7.1189836e+011 
Iter: 37 Ap: 1.69e-001 Pobj: -1.5153773e-001 Ad: 1.66e-001 Dobj: -7.1195303e+011 
Iter: 38 Ap: 2.11e-001 Pobj: -1.5226794e-001 Ad: 1.56e-001 Dobj: -7.1205845e+011 
Iter: 39 Ap: 2.21e-001 Pobj: -1.2269656e-001 Ad: 2.11e-001 Dobj: -7.1228465e+011 
Iter: 40 Ap: 2.26e-001 Pobj: -1.1381863e-001 Ad: 1.30e-001 Dobj: -7.1242018e+011 
Iter: 41 Ap: 1.90e-001 Pobj: -1.0404438e-001 Ad: 1.95e-001 Dobj: -7.1266116e+011 
Iter: 42 Ap: 2.09e-001 Pobj: -8.6792608e-002 Ad: 1.43e-001 Dobj: -7.1283776e+011 
Iter: 43 Ap: 2.56e-001 Pobj: -1.2611741e-001 Ad: 2.62e-001 Dobj: -7.1301279e+011
Iter: 44 Ap: 3.10e-001 Pobj: -1.5209232e-001 Ad: 2.77e-001 Dobj: -7.1300428e+011
Iter: 45 Ap: 2.26e-001 Pobj: -2.0416423e-001 Ad: 2.52e-001 Dobj: -7.1290687e+011
Iter: 46 Ap: 2.74e-001 Pobj: -2.9092366e-001 Ad: 3.54e-001 Dobj: -7.1275800e+011
Iter: 47 Ap: 9.36e-002 Pobj: -5.6264545e-001 Ad: 5.20e-002 Dobj: -7.1270029e+011
Iter: 48 Ap: 4.26e-002 Pobj: -6.6329424e-001 Ad: 4.91e-002 Dobj: -7.1259666e+011
Iter: 49 Ap: 3.63e-002 Pobj: -6.5939724e-001 Ad: 4.93e-002 Dobj: -7.1248753e+011
Declaring primal infeasibility.
Success: SDP is primal infeasible
Certificate of primal infeasibility: a'*y=-1.00000e+000, ||A'(y)-Z||=9.55843e-009
objective_value(model) = 0.6593972351693083
```
此时`L_max = 8`。

尝试一个之前一直想到但是没有来得及试的方法：指定初始值。然后CSDP还不支持初始值……

重新计算得到
```
CSDP 6.2.0
Iter:  0 Ap: 0.00e+000 Pobj:  0.0000000e+000 Ad: 0.00e+000 Dobj:  0.0000000e+000 
Iter:  1 Ap: 8.18e-002 Pobj: -1.3962147e-002 Ad: 9.43e-002 Dobj: -2.1809110e+010 
Iter:  2 Ap: 1.36e-001 Pobj: -4.2888011e-001 Ad: 6.88e-002 Dobj: -3.3221426e+010 
Iter:  3 Ap: 2.81e-002 Pobj: -1.0940146e+000 Ad: 2.17e-002 Dobj: -7.5508937e+010 
Iter:  4 Ap: 3.81e-002 Pobj: -8.7152731e-001 Ad: 7.41e-002 Dobj: -2.3816175e+011 
Iter:  5 Ap: 6.80e-002 Pobj:  1.4322870e+000 Ad: 5.05e-002 Dobj: -3.3485435e+011 
Iter:  6 Ap: 5.87e-002 Pobj:  3.5434138e+000 Ad: 6.03e-002 Dobj: -4.0310930e+011 
Iter:  7 Ap: 5.16e-002 Pobj:  3.6233577e+000 Ad: 4.89e-002 Dobj: -4.4682471e+011 
Iter:  8 Ap: 1.08e-001 Pobj:  4.1490790e+000 Ad: 9.53e-002 Dobj: -5.1735824e+011 
Iter:  9 Ap: 9.35e-002 Pobj:  2.4398432e+000 Ad: 8.56e-002 Dobj: -5.6546517e+011 
Iter: 10 Ap: 4.92e-002 Pobj:  2.4157495e+000 Ad: 5.79e-002 Dobj: -5.8655939e+011 
Iter: 11 Ap: 2.17e-001 Pobj:  1.3680341e+000 Ad: 1.11e-001 Dobj: -6.2108891e+011 
Iter: 12 Ap: 1.06e-001 Pobj:  1.6787347e+000 Ad: 1.37e-001 Dobj: -6.3750707e+011 
Iter: 13 Ap: 1.95e-001 Pobj:  1.5151120e+000 Ad: 9.04e-002 Dobj: -6.4720303e+011 
Iter: 14 Ap: 1.52e-001 Pobj:  3.3565862e-001 Ad: 1.63e-001 Dobj: -6.6165588e+011 
Iter: 15 Ap: 1.58e-001 Pobj:  7.6631334e-002 Ad: 1.44e-001 Dobj: -6.7205402e+011 
Iter: 16 Ap: 2.63e-001 Pobj: -3.0307305e-001 Ad: 1.90e-001 Dobj: -6.8199187e+011 
Iter: 17 Ap: 2.26e-001 Pobj: -8.9357966e-001 Ad: 1.75e-001 Dobj: -6.8839332e+011 
Iter: 18 Ap: 2.10e-001 Pobj: -6.7280867e-001 Ad: 1.11e-001 Dobj: -6.9183384e+011 
Iter: 19 Ap: 1.93e-001 Pobj: -8.8140377e-001 Ad: 1.66e-001 Dobj: -6.9566540e+011 
Iter: 20 Ap: 2.19e-001 Pobj: -8.0721903e-001 Ad: 7.31e-002 Dobj: -6.9702267e+011 
Iter: 21 Ap: 2.02e-001 Pobj: -6.8561861e-001 Ad: 2.40e-001 Dobj: -7.0163995e+011 
Iter: 22 Ap: 1.85e-001 Pobj: -7.2311332e-001 Ad: 1.74e-001 Dobj: -7.0378128e+011 
Iter: 23 Ap: 3.54e-001 Pobj: -7.1837629e-001 Ad: 2.75e-001 Dobj: -7.0714731e+011 
Iter: 24 Ap: 2.42e-001 Pobj: -3.7297316e-001 Ad: 2.09e-001 Dobj: -7.0850748e+011 
Iter: 25 Ap: 2.04e-001 Pobj: -2.6165012e-001 Ad: 1.95e-001 Dobj: -7.0913032e+011 
Iter: 26 Ap: 3.07e-001 Pobj: -1.6440458e-001 Ad: 1.94e-001 Dobj: -7.0999035e+011 
Iter: 27 Ap: 2.47e-001 Pobj: -1.4533222e-003 Ad: 1.88e-001 Dobj: -7.1029880e+011 
Iter: 28 Ap: 2.77e-001 Pobj:  2.5255903e-002 Ad: 2.05e-001 Dobj: -7.1068814e+011 
Iter: 29 Ap: 3.41e-001 Pobj: -9.4201483e-004 Ad: 2.37e-001 Dobj: -7.1083676e+011 
Iter: 30 Ap: 3.34e-001 Pobj: -2.0191262e-002 Ad: 2.04e-001 Dobj: -7.1095091e+011 
Iter: 31 Ap: 3.35e-001 Pobj:  6.0105422e-003 Ad: 2.64e-001 Dobj: -7.1100615e+011 
Iter: 32 Ap: 2.97e-001 Pobj: -1.7669119e-002 Ad: 2.46e-001 Dobj: -7.1107574e+011 
Iter: 33 Ap: 3.15e-001 Pobj: -1.2698698e-001 Ad: 1.93e-001 Dobj: -7.1120175e+011 
Iter: 34 Ap: 1.82e-001 Pobj: -2.1831267e-001 Ad: 9.58e-002 Dobj: -7.1125553e+011 
Iter: 35 Ap: 1.63e-001 Pobj: -2.0931588e-001 Ad: 1.03e-001 Dobj: -7.1131622e+011 
Iter: 36 Ap: 1.91e-001 Pobj: -1.7809097e-001 Ad: 1.25e-001 Dobj: -7.1136500e+011 
Iter: 37 Ap: 2.15e-001 Pobj: -1.4791290e-001 Ad: 2.14e-001 Dobj: -7.1147311e+011 
Iter: 38 Ap: 1.93e-001 Pobj: -1.4956594e-001 Ad: 1.60e-001 Dobj: -7.1160971e+011 
Iter: 39 Ap: 2.07e-001 Pobj: -1.0139010e-001 Ad: 1.99e-001 Dobj: -7.1189605e+011 
Iter: 40 Ap: 2.43e-001 Pobj: -1.0572125e-001 Ad: 1.36e-001 Dobj: -7.1204700e+011 
Iter: 41 Ap: 1.93e-001 Pobj: -1.1033933e-001 Ad: 1.85e-001 Dobj: -7.1232677e+011
Iter: 42 Ap: 2.34e-001 Pobj: -9.5756183e-002 Ad: 1.75e-001 Dobj: -7.1250214e+011
Iter: 43 Ap: 1.97e-001 Pobj: -1.1980201e-001 Ad: 2.71e-001 Dobj: -7.1263298e+011
Iter: 44 Ap: 2.42e-001 Pobj: -1.3890931e-001 Ad: 2.34e-001 Dobj: -7.1264511e+011
Iter: 45 Ap: 2.74e-001 Pobj: -1.8944747e-001 Ad: 3.09e-001 Dobj: -7.1260420e+011
Iter: 46 Ap: 1.89e-001 Pobj: -2.6079582e-001 Ad: 2.42e-001 Dobj: -7.1244221e+011
Iter: 47 Ap: 1.58e-001 Pobj: -3.5200318e-001 Ad: 1.93e-001 Dobj: -7.1229941e+011
Declaring primal infeasibility.
Success: SDP is primal infeasible
Certificate of primal infeasibility: a'*y=-1.00000e+000, ||A'(y)-Z||=9.88450e-009
objective_value(model) = 0.3520031810148794
```
可以看到从中间某一步出发两次计算就不一致了。

idea：解析计算出所有的关联函数，然后验证它代入此处的问题是否是feasible的。如果是，那么我的程序的毛病就是纯粹技术性的，如果不是，我肯定就哪里写错了。

## 2022.3.19

先安装了`COSMO.jl`，更换一下求解器试试。

解析计算了一些结果(`calculate-all-correlation-functions-in-xp-oscillator-for-benchmark-1.nb`)
得到
```
{{xpOpString[p, p] -> 0.808702 + 0. I, xpOpString[x, p, p] -> 0, 
  xpOpString[x, x, p, p] -> -0.186501 + 0. I, 
  xpOpString[x, x, x, p, p] -> 0, 
  xpOpString[x, x, x, x, p, p] -> -0.595472 + 0. I, 
  xpOpString[x, x, x, x, x, p, p] -> 0, 
  xpOpString[x, x, x, x, x, x, p, p] -> -1.44281 + 0. I}}
```

这次倒是没有不可解了，又不收敛了。
```
------------------------------------------------------------------
          COSMO v0.8.2 - A Quadratic Objective Conic Solver
                         Michael Garstka
                University of Oxford, 2017 - 2021
------------------------------------------------------------------

Problem:  x ∈ R^{13779},
          constraints: A ∈ R^{26614x13779} (82577 nnz),
          matrix size to factor: 40393x40393,
          Floating-point precision: Float64
Sets:     ZeroSe of dim: 13411
          DensePsdConeTriangl of dim: 13203 (162x162)
Settings: ϵ_abs = 1.0e-05, ϵ_rel = 1.0e-05,
          ϵ_prim_inf = 1.0e-04, ϵ_dual_inf = 1.0e-04,
          ρ = 0.1, σ = 1e-06, α = 1.6,
          max_iter = 5000,
          scaling iter = 10 (on),
          check termination every 25 iter,
          check infeasibility every 40 iter,
          KKT system solver: QDLDL
Acc:      Anderson Type2{QRDecomp},
          Memory size = 15, RestartedMemory,
          Safeguarded: true, tol: 2.0
Setup Time: 814.6ms

Iter:   Objective:      Primal Res:     Dual Res:       Rho:
1       -1.6145e-01     2.6544e+04      1.0865e+02      1.0000e-01
25       8.6660e-02     3.9350e+03      1.0342e+02      1.0000e-01
50      -1.7859e-02     4.2178e+03      4.8657e+01      1.0000e-01
75       2.2133e-02     3.2869e+03      4.2510e+01      1.0000e-01
100      5.2120e-02     1.9585e+03      1.5275e+01      1.6514e-02
125      1.0539e-02     3.6843e+03      1.6217e+00      2.0670e-03
150     -1.9105e-01     7.8851e+03      1.2819e+01      2.0670e-03
175     -1.5812e-01     1.0420e+04      1.1964e+01      2.0670e-03
200     -2.8632e-01     3.3036e+04      1.0004e+01      2.0670e-03
225     -1.0203e+00     2.2032e+05      2.7876e+01      2.0610e-04
250     -2.5277e+00     3.5177e+05      1.3187e+01      1.1020e-05
275     -2.7546e+01     2.9427e+06      5.7900e+01      1.1020e-05
300      1.4341e+01     1.6311e+07      1.1605e+02      1.0000e-06
325     -6.7725e+01     1.6188e+07      6.7111e+01      1.0000e-06
350     -1.4972e+02     1.7605e+07      1.1260e+02      1.0000e-06
375     -1.4214e+02     1.4745e+07      1.4506e+02      1.0000e-06
400     -1.0606e+02     1.7220e+07      1.4594e+02      1.0000e-06
425     -8.5678e+01     1.4931e+07      1.4085e+02      1.0000e-06
450     -1.1112e+02     1.4907e+07      1.4209e+02      1.0000e-06
475     -1.4832e+02     1.4131e+07      1.5083e+02      1.0000e-06
500     -1.7313e+02     1.4663e+07      1.4471e+02      1.0000e-06
525     -2.0047e+02     1.5522e+07      1.4642e+02      1.0000e-06
550     -1.8126e+02     9.6146e+06      1.4481e+02      1.0000e-06
575     -2.0708e+02     1.1391e+07      1.2719e+02      1.0000e-06
600     -2.3476e+02     9.8822e+06      8.4573e+01      1.0000e-06
625     -2.6161e+02     7.6195e+06      7.5813e+01      1.0000e-06
650     -2.7654e+02     1.0188e+07      5.1954e+01      1.0000e-06
675     -2.6234e+02     6.1881e+06      5.7449e+01      1.0000e-06
700     -2.4793e+02     5.1644e+06      6.2402e+01      1.0000e-06
725     -2.3552e+02     1.1234e+07      6.0566e+01      1.0000e-06
750     -2.2169e+02     6.8768e+06      4.5006e+01      1.0000e-06
775     -2.2912e+02     5.6048e+06      4.5746e+01      1.0000e-06
800     -2.2619e+02     1.1086e+07      5.2272e+01      1.0000e-06
825     -2.2046e+02     6.4686e+06      5.9063e+01      1.0000e-06
850     -2.0967e+02     3.3002e+06      6.3651e+01      1.0000e-06
875     -2.0629e+02     2.9677e+06      7.3216e+01      1.0000e-06
900     -1.9123e+02     6.9276e+06      8.9149e+01      1.0000e-06
925     -1.9299e+02     3.2376e+06      8.8018e+01      1.0000e-06
950     -1.9876e+02     4.5456e+06      8.9052e+01      1.0000e-06
975     -1.9984e+02     6.5809e+06      8.9809e+01      1.0000e-06
1000    -2.0036e+02     3.4840e+06      8.9282e+01      1.0000e-06
1025    -1.8941e+02     5.7877e+06      9.3482e+01      1.0000e-06
1050    -1.7552e+02     8.2936e+06      9.3508e+01      1.0000e-06
1075    -1.7418e+02     4.6710e+06      9.2227e+01      1.0000e-06
1100    -1.7109e+02     2.5637e+06      9.0228e+01      1.0000e-06
1125    -1.6926e+02     7.8451e+06      8.5086e+01      1.0000e-06
1150    -1.7048e+02     5.0177e+06      8.4248e+01      1.0000e-06
1175    -1.6165e+02     2.3359e+06      8.4802e+01      1.0000e-06
1200    -1.3658e+02     7.1121e+06      8.1617e+01      1.0000e-06
1225    -1.2807e+02     5.3249e+06      7.8584e+01      1.0000e-06
1250    -1.2398e+02     3.1177e+06      7.6562e+01      1.0000e-06
1275    -1.2232e+02     3.3343e+06      7.1100e+01      1.0000e-06
1300    -1.1995e+02     3.6417e+06      6.8756e+01      1.0000e-06
1325    -1.1238e+02     1.5598e+06      6.6756e+01      1.0000e-06
1350    -9.5565e+01     3.5903e+06      6.5733e+01      1.0000e-06
1375    -8.5177e+01     2.9803e+06      6.8096e+01      1.0000e-06
1400    -8.0302e+01     1.3060e+06      6.7344e+01      1.0000e-06
1425    -7.9799e+01     3.1009e+06      6.4733e+01      1.0000e-06
1450    -7.9414e+01     3.7008e+06      6.1095e+01      1.0000e-06
1475    -7.6824e+01     2.3414e+06      5.9825e+01      1.0000e-06
1500    -7.0628e+01     9.1390e+05      5.8037e+01      1.0000e-06
1525    -5.1073e+01     4.8671e+06      5.5741e+01      1.0000e-06
1550    -4.6168e+01     3.9369e+06      5.4957e+01      1.0000e-06
1575    -4.4896e+01     3.1309e+06      5.4030e+01      1.0000e-06
1600    -4.4247e+01     1.1885e+06      5.2120e+01      1.0000e-06
1625    -4.7478e+01     3.7319e+06      4.6059e+01      1.0000e-06
1650    -4.7642e+01     3.3940e+06      4.2055e+01      1.0000e-06
1675    -4.2731e+01     8.0081e+05      3.9079e+01      1.0000e-06
1700    -3.2865e+01     1.6932e+06      3.4754e+01      1.0000e-06
1725    -2.8022e+01     1.5382e+06      3.0257e+01      1.0000e-06
1750    -2.5227e+01     7.1300e+05      2.7729e+01      1.0000e-06
1775    -2.5516e+01     1.2635e+06      2.2620e+01      1.0000e-06
1800    -2.5306e+01     1.7188e+06      2.4323e+01      1.0000e-06
1825    -2.3049e+01     9.4920e+05      2.4577e+01      1.0000e-06
1850    -1.8155e+01     1.5760e+06      2.5442e+01      1.0000e-06
1875    -1.3003e+01     2.1789e+06      2.5894e+01      1.0000e-06
1900    -1.1385e+01     1.4193e+06      2.6229e+01      1.0000e-06
1925    -1.2252e+01     5.7776e+05      2.7228e+01      1.0000e-06
1950    -1.4887e+01     2.0467e+06      3.1700e+01      1.0000e-06
1975    -1.5242e+01     1.3906e+06      3.3705e+01      1.0000e-06
2000    -1.4436e+01     7.0826e+05      3.5231e+01      1.0000e-06
2025    -8.9350e+00     2.7505e+06      3.8606e+01      1.0000e-06
2050    -8.3847e+00     2.2328e+06      4.0102e+01      1.0000e-06
2075    -8.3508e+00     2.0323e+06      4.0490e+01      1.0000e-06
2100    -9.5369e+00     1.1128e+06      4.1816e+01      1.0000e-06
2125    -1.2813e+01     1.1103e+06      4.5113e+01      1.0000e-06
2150    -1.7442e+01     2.2868e+06      4.7365e+01      1.0000e-06
2175    -1.7139e+01     1.1933e+06      4.6948e+01      1.0000e-06
2200    -1.6578e+01     1.0305e+06      4.7407e+01      1.0000e-06
2225    -1.5409e+01     1.0793e+06      4.7650e+01      1.0000e-06
2250    -1.6859e+01     5.3375e+05      4.8622e+01      1.0000e-06
2275    -2.2051e+01     1.9853e+06      5.0652e+01      1.0000e-06
2300    -2.3780e+01     2.0210e+06      5.0580e+01      1.0000e-06
2325    -2.3777e+01     1.4735e+06      5.0084e+01      1.0000e-06
2350    -2.5058e+01     1.1969e+06      4.9862e+01      1.0000e-06
2375    -2.4265e+01     1.5915e+06      4.7825e+01      1.0000e-06
2400    -2.3015e+01     1.9499e+06      4.7083e+01      1.0000e-06
2425    -2.3312e+01     1.3889e+06      4.6485e+01      1.0000e-06
2450    -2.5204e+01     6.0393e+05      4.6802e+01      1.0000e-06
2475    -3.1212e+01     2.7515e+06      4.5885e+01      1.0000e-06
2500    -3.2011e+01     2.1212e+06      4.4936e+01      1.0000e-06
2525    -3.1503e+01     1.8462e+06      4.4539e+01      1.0000e-06
2550    -3.1964e+01     6.6416e+05      4.2338e+01      1.0000e-06
2575    -3.0747e+01     2.0771e+06      3.9297e+01      1.0000e-06
2600    -3.0233e+01     1.7090e+06      3.8227e+01      1.0000e-06
2625    -3.0639e+01     9.0436e+05      3.7761e+01      1.0000e-06
2650    -3.4892e+01     1.7814e+06      3.5902e+01      1.0000e-06
2675    -3.6664e+01     1.7468e+06      3.4314e+01      1.0000e-06
2700    -3.6357e+01     1.1137e+06      3.2458e+01      1.0000e-06
2725    -3.6511e+01     7.4746e+05      3.1686e+01      1.0000e-06
2750    -3.4355e+01     1.2695e+06      2.8278e+01      1.0000e-06
2775    -3.4784e+01     7.5102e+05      2.8671e+01      1.0000e-06
2800    -3.6113e+01     1.1667e+06      2.4648e+01      1.0000e-06
2825    -3.8642e+01     1.8877e+06      2.2457e+01      1.0000e-06
2850    -3.8157e+01     1.0026e+06      2.1035e+01      1.0000e-06
2875    -3.7561e+01     4.3785e+05      2.0655e+01      1.0000e-06
2900    -3.6588e+01     5.1822e+05      1.9018e+01      1.0000e-06
2925    -3.6082e+01     3.6263e+05      1.6546e+01      1.0000e-06
2950    -3.5773e+01     6.3119e+05      1.6678e+01      1.0000e-06
2975    -3.5834e+01     4.2474e+05      1.6948e+01      1.0000e-06
3000    -3.5203e+01     5.4711e+05      1.7706e+01      1.0000e-06
3025    -3.5069e+01     3.5481e+05      1.8294e+01      1.0000e-06
3050    -3.4636e+01     3.6518e+05      1.8649e+01      1.0000e-06
3075    -3.3528e+01     3.6893e+05      1.9071e+01      1.0000e-06
3100    -3.2086e+01     4.5241e+05      1.9451e+01      1.0000e-06
3125    -3.1690e+01     3.8809e+05      1.9693e+01      1.0000e-06
3150    -3.0908e+01     3.7164e+05      1.9812e+01      1.0000e-06
3175    -3.0323e+01     7.1646e+05      1.9952e+01      1.0000e-06
3200    -2.9875e+01     4.5787e+05      1.9967e+01      1.0000e-06
3225    -2.8671e+01     3.9602e+05      1.9388e+01      1.0000e-06
3250    -2.6007e+01     9.1586e+05      1.8795e+01      1.0000e-06
3275    -2.6229e+01     5.1864e+05      1.8871e+01      1.0000e-06
3300    -2.5483e+01     4.0072e+05      1.8708e+01      1.0000e-06
3325    -2.5394e+01     3.7429e+05      1.8321e+01      1.0000e-06
3350    -2.5339e+01     7.2098e+05      1.7836e+01      1.0000e-06
3375    -2.4106e+01     4.3320e+05      1.7479e+01      1.0000e-06
3400    -2.2154e+01     4.4845e+05      1.6703e+01      1.0000e-06
3425    -2.0943e+01     1.0429e+06      1.5315e+01      1.0000e-06
3450    -2.0711e+01     6.7017e+05      1.5179e+01      1.0000e-06
3475    -2.0827e+01     4.9231e+05      1.4910e+01      1.0000e-06
3500    -2.0791e+01     9.5022e+05      1.4109e+01      1.0000e-06
3525    -2.1386e+01     1.0681e+06      1.3746e+01      1.0000e-06
3550    -2.1292e+01     8.0279e+05      1.3716e+01      1.0000e-06
3575    -2.0423e+01     3.9166e+05      1.3970e+01      1.0000e-06
3600    -1.9070e+01     7.6172e+05      1.4878e+01      1.0000e-06
3625    -1.9205e+01     6.1196e+05      1.5461e+01      1.0000e-06
3650    -1.8574e+01     4.2474e+05      1.5683e+01      1.0000e-06
3675    -1.9601e+01     7.0766e+05      1.6607e+01      1.0000e-06
3700    -2.0230e+01     1.0341e+06      1.7290e+01      1.0000e-06
3725    -2.0979e+01     7.5016e+05      1.7336e+01      1.0000e-06
3750    -2.0257e+01     4.9604e+05      1.7534e+01      1.0000e-06
3775    -1.8299e+01     1.1492e+06      1.8173e+01      1.0000e-06
3800    -1.8311e+01     1.1624e+06      1.8314e+01      1.0000e-06
3825    -1.8511e+01     7.3889e+05      1.8417e+01      1.0000e-06
3850    -1.9336e+01     3.6616e+05      1.8853e+01      1.0000e-06
3875    -2.0298e+01     5.4809e+05      1.9177e+01      1.0000e-06
3900    -2.1779e+01     3.8413e+05      1.9427e+01      1.0000e-06
3925    -2.1474e+01     3.7393e+05      1.9242e+01      1.0000e-06
3950    -2.1313e+01     5.4570e+05      1.9340e+01      1.0000e-06
3975    -2.0689e+01     5.2317e+05      1.9439e+01      1.0000e-06
4000    -2.1679e+01     3.8009e+05      1.9542e+01      1.0000e-06
4025    -2.3126e+01     7.6652e+05      1.9857e+01      1.0000e-06
4050    -2.3710e+01     6.9333e+05      1.9750e+01      1.0000e-06
4075    -2.4454e+01     4.2314e+05      1.9610e+01      1.0000e-06
4100    -2.4276e+01     3.9808e+05      1.9081e+01      1.0000e-06
4125    -2.4795e+01     1.3214e+06      1.8219e+01      1.0000e-06
4150    -2.4546e+01     1.1746e+06      1.8394e+01      1.0000e-06
4175    -2.4469e+01     1.1530e+06      1.8679e+01      1.0000e-06
4200    -2.4224e+01     1.1188e+06      1.8596e+01      1.0000e-06
4225    -2.4499e+01     1.0241e+06      1.8267e+01      1.0000e-06
4250    -2.4497e+01     8.4866e+05      1.8463e+01      1.0000e-06
4275    -2.6403e+01     3.7473e+05      1.8436e+01      1.0000e-06
4300    -2.7256e+01     4.9165e+05      1.8409e+01      1.0000e-06
4325    -2.8173e+01     4.7511e+05      1.8031e+01      1.0000e-06
4350    -2.8464e+01     3.8183e+05      1.7730e+01      1.0000e-06
4375    -2.8718e+01     5.0224e+05      1.6594e+01      1.0000e-06
4400    -2.8667e+01     6.6969e+05      1.6217e+01      1.0000e-06
4425    -2.8957e+01     4.6308e+05      1.6095e+01      1.0000e-06
4450    -3.0012e+01     3.8532e+05      1.5890e+01      1.0000e-06
4475    -3.1978e+01     9.2135e+05      1.5170e+01      1.0000e-06
4500    -3.2431e+01     7.1787e+05      1.4740e+01      1.0000e-06
4525    -3.2007e+01     4.3350e+05      1.4002e+01      1.0000e-06
4550    -3.1689e+01     7.0903e+05      1.2926e+01      1.0000e-06
4575    -3.1489e+01     1.0245e+06      1.2316e+01      1.0000e-06
4600    -3.2014e+01     9.2615e+05      1.2120e+01      1.0000e-06
4625    -3.1749e+01     7.7482e+05      1.2142e+01      1.0000e-06
4650    -3.2979e+01     4.6520e+05      1.1904e+01      1.0000e-06
4675    -3.2517e+01     3.9758e+05      1.1854e+01      1.0000e-06
4700    -3.4099e+01     5.2907e+05      1.1066e+01      1.0000e-06
4725    -3.4440e+01     6.1480e+05      1.0712e+01      1.0000e-06
4750    -3.4330e+01     4.3538e+05      1.0165e+01      1.0000e-06
4775    -3.4699e+01     4.2150e+05      9.4425e+00      1.0000e-06
4800    -3.3202e+01     7.5848e+05      1.4530e+01      1.0000e-06
4825    -3.4034e+01     6.2543e+05      8.0730e+00      1.0000e-06
4850    -3.4218e+01     4.6986e+05      8.0888e+00      1.0000e-06
4875    -3.4629e+01     4.0084e+05      7.8625e+00      1.0000e-06
4900    -3.5727e+01     6.8322e+05      6.5542e+00      1.0000e-06
4925    -3.5858e+01     4.5778e+05      6.5247e+00      1.0000e-06
4950    -3.5467e+01     4.0434e+05      6.5011e+00      1.0000e-06
4975    -3.4185e+01     9.5330e+05      8.4012e+00      1.0000e-06
5000    -3.4141e+01     7.8652e+05      6.8870e+00      1.0000e-06

------------------------------------------------------------------
>>> Results
Status: Max_iter_reached
Iterations: 5000
Optimal objective: -34.14
Runtime: 99.763s (99763.0ms)
```

那好吧，我再试试`SCS`。

看起来也不是太行的样子。放弃了……

将`L_max`降低到5，观察SCS会给出什么结果。

结果如下：
```
----------------------------------------------------------------------------
Status: Infeasible/Inaccurate
Hit max_iters, solution may be inaccurate, returning best found solution.
Timing: Solve time: 1.54e+002s
        Lin-sys: avg # CG iterations: 453.00, avg solve time: 2.33e-002s
        Cones: avg projection time: 6.19e-003s
        Acceleration: avg step time: 1.10e-003s
----------------------------------------------------------------------------
Certificate of primal infeasibility:
dist(y, K*) = 1.4565e-009
|A'y|_2 * |b|_2 = 1.8505e+000
b'y = -1.0000
============================================================================
```
再换用COSMO，得到
```
------------------------------------------------------------------
          COSMO v0.8.2 - A Quadratic Objective Conic Solver
                         Michael Garstka
                University of Oxford, 2017 - 2021
------------------------------------------------------------------

Problem:  x ∈ R^{2868},
          constraints: A ∈ R^{5311x2868} (13134 nnz),
          matrix size to factor: 8179x8179,
          Floating-point precision: Float64
Sets:     ZeroSe of dim: 2683
          DensePsdConeTriangl of dim: 2628 (72x72)
Settings: ϵ_abs = 1.0e-05, ϵ_rel = 1.0e-05,
          ϵ_prim_inf = 1.0e-04, ϵ_dual_inf = 1.0e-04,
          ρ = 0.1, σ = 1e-06, α = 1.6,
          max_iter = 5000,
          scaling iter = 10 (on),
          check termination every 25 iter,
          check infeasibility every 40 iter,
          KKT system solver: QDLDL
Acc:      Anderson Type2{QRDecomp},
          Memory size = 15, RestartedMemory,
          Safeguarded: true, tol: 2.0
Setup Time: 698.88ms

Iter:   Objective:      Primal Res:     Dual Res:       Rho:
1       -3.9336e-01     6.9832e+01      8.0683e-01      1.0000e-01
25       1.5816e-01     2.8502e+00      1.7009e-01      1.0000e-01
50       1.7040e-01     2.7070e+00      1.8386e-01      1.2427e-02
75       3.3340e-01     3.5199e+00      1.4083e-01      1.2427e-02

------------------------------------------------------------------
>>> Results
Status: Primal_infeasible
Iterations: 87 (incl. 6 safeguarding iter)
Optimal objective: Inf
Runtime: 3.947s (3947.0ms)
```
因此程序多半实在是写错了。

唉，现在只能老实检查了。

建立`calculate-all-correlation-functions-in-xp-oscillator-for-benchmark-2.nb`，用于系统性地计算L_max=5的情况。
其实还有一件事可以做，就是在Mathematica里面画图。不知道存在量词这一套能不能用在`xpOpString[x ,x ,p]`这样的东西上面。

下面可以做的事情：
- 检查约束情况
  - 计算`L_max = 5`的Mathematica benchmark；幂次最高为$x^{10} p^{10}$
  - 但是这里就有一个问题，就是递归计算出来，高次关联函数都是很大的，是不是会不准确等等……
  - 不过还是先检查约束是否可以满足吧
- 用Mathematica画图
- 说起来我感觉我的思维方式被限制住了，为什么不能直接用解薛定谔方程得到参考数值？

## 2022.3.20

在`stationary-schrodinger-2.jl`中计算$\lang x^m p^n \rang$。

汇总目前通过严格做numerical bootstrap得到的数据：
[此处](#2022314)的$\lang x^{2n} \rang$
```
0.301138, 0.253782, 0.343358, 0.598177, 1.31284, 3.40689, 9.75835, 
32.4658, 117.962, 451.727
```
[此处](#2022319)的
```
{{xpOpString[p, p] -> 0.808702 + 0. I, xpOpString[x, p, p] -> 0, 
  xpOpString[x, x, p, p] -> -0.186501 + 0. I, 
  xpOpString[x, x, x, p, p] -> 0, 
  xpOpString[x, x, x, x, p, p] -> -0.595472 + 0. I, 
  xpOpString[x, x, x, x, x, p, p] -> 0, 
  xpOpString[x, x, x, x, x, x, p, p] -> -1.44281 + 0. I}}
```
从`calculate-all-correlation-functions-in-xp-oscillator-for-benchmark-1.nb`读出一些之前没有记录的数据：
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

用它们和`stationary-schrodinger-2.jl`的输出相对比：$\lang x^n \rang$
```
     1.0000000000000004 + 0.0im
 -9.672537983335926e-15 + 0.0im
    0.30579492589833185 + 0.0im
 -7.419859045292633e-15 + 0.0im
     0.2602109595243272 + 0.0im
 -9.215911011408188e-15 + 0.0im
    0.34723077229789145 + 0.0im
  -1.52123002337647e-14 + 0.0im
     0.6163426616900493 + 0.0im
 -3.171550048562737e-14 + 0.0im
     1.3459312154711665 + 0.0im
```
$\lang x^n p \rang$
```
0.0 + 0.008262813010199016im
 0.0 + 0.4999173718698987im
 0.0 - 0.0018168737481579495im
 0.0 + 0.45864691158556825im
 0.0 - 0.006012850492717911im
 0.0 + 0.6505222900588602im
 0.0 - 0.014788289221550878im
 0.0 + 1.2154322376908027im
 0.0 - 0.03943614114283746im
 0.0 + 2.774174630029892im
 0.0 - 0.11712209240572102im
```
$\lang x^n p^2 \rang$
```
 0.8258946223615176 + 0.0im
  -0.024780705457444523 + 0.0im
   -0.18139193503898404 + 0.0im
  -0.013641080401041665 + 0.0im
     -0.600674817668086 + 0.0im
 -0.0015412857505892742 + 0.0im
    -1.4777566743789248 + 0.0im
    0.03732157515563892 + 0.0im
     -3.941977898035117 + 0.0im
    0.18966713392589257 + 0.0im
     -11.71143017592804 + 0.0im
```
大致对得上，虽然误差很有些大。

ok，现在数据有了，可以试着做测试了。输出数据保存在`xpopstr_expected_ode-dx-0.02-l-10.jld2`。

在`jump-oscillator-2-benchmark-with-ode.jl`中讨论我自己写的constraint是否是正确的，使用`xpopstr_expected_ode-dx-0.02-l-10.jld2`中的数据做benchmark。

首先尝试直接将优化variable替换成真实值，看看约束是不是成立；可是这第一步就出了奇怪的问题。运行如下代码：
```julia
fix(xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :real), 0.298)
fix(xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :imag), 0.0)
fix(xpopstr_expected_real_imag_parts(xpopstr_index(4, 0), :real), 0.258)
fix(xpopstr_expected_real_imag_parts(xpopstr_index(4, 0), :imag), 0.0)

optimize!(model)

value(xpopstr_expected_real_imag_parts(xpopstr_index(4, 0), :real))
```
发现输出结果是`0.21942580250281907`。也就是说拿`fix`根本fix不了变量。有毒。

退而求其次，观察约束的系数是否正确。

采用如下代码构建测试数据：
```julia
using JLD2
working_path = "D:\\Projects\\numerical-boostrap\\oscillator-simple-prototype\\"
@load working_path * "xpopstr_expected_ode-dx-0.02-l-10.jld2" xpopstr_expected_ode 

xpopstr_expected_ope_value = zeros(2xpopspace_dim)

for xpopstr_idx in 1 : xpopspace_dim
    x_power = index_to_xpower(xpopstr_idx)
    p_power = index_to_ppower(xpopstr_idx)
    expected_value = xpopstr_expected_ode[x_power, p_power]
    xpopstr_expected_ope_value[2xpopstr_idx - 1] = real(expected_value)
    xpopstr_expected_ope_value[2xpopstr_idx] = imag(expected_value)
end

xpopstr_expected_real_imag_parts_ope_value(i, real_or_imag) = begin
    if real_or_imag == :real
        return xpopstr_expected_ope_value[2i - 1]
    end
    if real_or_imag == :imag
        return xpopstr_expected_ope_value[2i]
    end
end

variable_list_real_ope_value = 
    [xpopstr_expected_real_imag_parts_ope_value(i, :real) * I22 for i in 1 : xpopspace_dim]
variable_list_imag_ope_value = 
    [xpopstr_expected_real_imag_parts_ope_value(i, :imag) * Im22 for i in 1 : xpopspace_dim]
xpopstr_basis_real_ope_value = OffsetArray([I22, variable_list_real_ope_value...], xpopspace_index_range)
xpopstr_basis_imag_ope_value = OffsetArray([Im22, variable_list_imag_ope_value...], xpopspace_index_range)

function complex_to_mat_ope_value(coefficients)
    real_part = transpose(real(coefficients))
    imag_part = transpose(imag(coefficients))
    real_part_mat_version = map(x -> x * I22, real_part)
    imag_part_mat_version = map(x -> x * Im22, imag_part)
    (real_part_mat_version + imag_part_mat_version) * (xpopstr_basis_real_ope_value + xpopstr_basis_imag_ope_value)
end

complex_to_mat_ope_value(xpopstr_xp_power(4, 0))
```
测试$\lang [H, O] \rang = 0$的诊断代码：
```julia
for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = real(cons)
    cons_imag = imag(cons)
    lhs = complex_to_mat(cons)
    if cons_real == cons_imag == zero_xpopstr
        continue
    end
    if cons_real == zero_xpopstr
        println(lhs[1, 2])
    elseif cons_imag == zero_xpopstr
        println(lhs[1, 1])
    else
        println(lhs[1, 1])
        println(lhs[1, 2])
    end
end
```
输出为
```
-4.902451214784239e-14
-0.0006943306158411211
2.0
0.0
-0.32820777442942706
0.0
0.0006444451709375354
-0.17131613487814604
0.0
0.0
0.020025042775518642
```
发现三处可疑的地方。改用如下诊断代码，定位可疑之处：
```julia
for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = real(cons)
    cons_imag = imag(cons)
    lhs = complex_to_mat_ope_value(cons)
    if cons_real == cons_imag == zero_xpopstr
        continue
    end
    if cons_real == zero_xpopstr
        println("$x_power $p_power   $(lhs[1, 2])")
    elseif cons_imag == zero_xpopstr
        println("$x_power $p_power   $(lhs[1, 1])")
    else
        println("$x_power $p_power   $(lhs[1, 2])")
        println("$x_power $p_power   $(lhs[1, 1])")
    end
end
```
得到
```
0 1   -4.902451214784239e-14
0 2   2.0
0 2   -0.0006943306158411211
0 3   -0.32820777442942706  
0 3   0.0
1 0   0.0
1 1   0.0006444451709375354
1 2   0.0
1 2   -0.17131613487814604
1 3   0.020025042775518642
1 3   0.0
```
离谱之处：$p^2$, $p^3$和$x p$。

以$p^2$为例。运行
```julia
xpopstr_stringify(comm_with_ham(xpopstr_xp_power(0, 2)))
```
得到
```
- 2.0 - 4.0im x p - 12.0 x^2 - 8.0im x^3 p
```
与手动计算一致。将所有变量用另行计算的值代替，即
```julia
- 2.0 - 4.0im * (0.5im) - 12.0 * 0.301138 - 8.0im * (0.451707im)

##

- 2.0 - 4.0im * (0.4999173718698987im) - 12.0 * (0.30579492589833185) - 8.0im * (0.45864691158556825im)
```
得到
```
0.0 - 0.0im

-0.0006943306158411211 - 0.0im
```
和前述自动计算的明显不一样啊。另一方面，
```julia
complex_to_mat_ope_value(comm_with_ham(xpopstr_xp_power(0, 2)))
```
给出
```
2×2 Array{Float64,2}:
 -0.000694331   2.0
 -2.0          -0.000694331
```

在下面的代码中出现了惊人一幕：
```julia
coefficients = comm_with_ham(xpopstr_xp_power(0, 2))
real_part = transpose(real(coefficients))
imag_part = transpose(imag(coefficients))
real_part_mat_version = map(x -> x * I22, real_part)
imag_part_mat_version = map(x -> x * Im22, imag_part)
(real_part_mat_version + imag_part_mat_version) * (xpopstr_basis_real_ope_value + xpopstr_basis_imag_ope_value)
```
运行结果为
```
julia> coefficients[xpopstr_index(2, 0)]
-12.0 + 0.0im

julia> real_part[xpopstr_index(2, 0)]
0.0
```
我好像有点知道了，julia里面OffsetArray的转置是有问题的。考虑以下情况：
```
julia> transpose(coefficients)[xpopstr_index(2, 0)]
0.0 + 0.0im

julia> transpose(coefficients)[1, xpopstr_index(2, 0)]
-12.0 + 0.0im
```
更加清楚的是如下代码：
```
julia> transpose(OffsetArray([1, 2, 3], 0:2))[1]
1

julia> transpose(OffsetArray([1, 2, 3], 0:2))[2]
2

julia> transpose(OffsetArray([1, 2, 3], 0:2))[3]
3
```
总之有下面的结果：
```
julia> xpopstr_basis_real_ope_value[xpopstr_index(2, 0)]
2×2 Array{Float64,2}:
 0.305795  0.0
 0.0       0.305795

julia> xpopstr_basis_imag_ope_value[xpopstr_index(2, 0)]
2×2 Array{Float64,2}:
 0.0  -0.0
 0.0   0.0

julia> coefficients[xpopstr_index(2, 0)]
-12.0 + 0.0im

julia> real_part[xpopstr_index(2, 0)]
0.0
```

现在的问题是，这是否会影响乘法计算。为此考虑如下诊断代码：
```julia
coefficients = xpopstr_xp_power(2, 0) 
real_part = transpose(real(coefficients))
imag_part = transpose(imag(coefficients))
real_part_mat_version = map(x -> x * I22, real_part)
imag_part_mat_version = map(x -> x * Im22, imag_part)
(real_part_mat_version + imag_part_mat_version) * (xpopstr_basis_real_ope_value + xpopstr_basis_imag_ope_value)
```
这给出
```
2×2 Array{Float64,2}:
 0.305795  0.0
 0.0       0.305795
```
这个时候就没有问题了……

我现在感觉，问题并不在这里的乘法上面。至少有一处地方是有错的，就是
```julia
xpopstr_basis_imag = OffsetArray([Im22, variable_list_imag...], xpopspace_index_range)
```
应该改成
```julia
xpopstr_basis_imag = OffsetArray([O22, variable_list_imag...], xpopspace_index_range)
```

解决此问题后，得到
```
Iter: 60 Ap: 1.00e+000 Pobj: -2.3212105e+000 Ad: 1.00e+000 Dobj: -2.7433469e+000 
Iter: 61 Ap: 1.00e+000 Pobj: -2.3214312e+000 Ad: 1.00e+000 Dobj: -2.7436446e+000 
Iter: 62 Ap: 1.00e+000 Pobj: -2.3214572e+000 Ad: 1.00e+000 Dobj: -2.7437136e+000 
Lack of progress.  Giving up!
Partial Success: SDP solved with reduced accuracy
Primal objective value: -2.7259246e+000
Dual objective value: -3.7969216e+000
Relative primal infeasibility: 1.89e-006
Relative dual infeasibility: 7.10e-007
Real Relative Gap: -1.42e-001
XZ Relative Gap: 3.53e-006
DIMACS error measures: 2.72e-006 0.00e+000 1.58e-006 0.00e+000 -1.42e-001 3.53e-006
objective_value(model) = 2.7259246483686406
2.7259246483686406
```
增大`L_max`以后并没有什么改善；实际上事情变得更糟糕了：
```
CSDP 6.2.0
Iter:  0 Ap: 0.00e+000 Pobj:  0.0000000e+000 Ad: 0.00e+000 Dobj:  0.0000000e+000 
Iter:  1 Ap: 8.80e-001 Pobj: -3.8497934e-001 Ad: 8.24e-001 Dobj: -2.3755443e+005 
Iter:  2 Ap: 7.66e-001 Pobj: -1.1723146e-001 Ad: 7.50e-001 Dobj: -8.7864365e+005 
Iter:  3 Ap: 4.74e-001 Pobj: -1.3935154e-001 Ad: 5.62e-001 Dobj: -4.5625669e+006 
Iter:  4 Ap: 4.50e-001 Pobj: -1.6926176e-001 Ad: 4.56e-001 Dobj: -7.3760066e+006 
Iter:  5 Ap: 3.44e-001 Pobj: -2.1358992e-001 Ad: 2.97e-001 Dobj: -8.4295899e+006 
Iter:  6 Ap: 1.53e-001 Pobj: -2.2305444e-001 Ad: 1.42e-001 Dobj: -1.3591702e+007 
Iter:  7 Ap: 2.62e-001 Pobj: -2.6397489e-001 Ad: 2.43e-001 Dobj: -1.9624439e+007 
Iter:  8 Ap: 1.17e-001 Pobj: -2.4948540e-001 Ad: 1.34e-001 Dobj: -2.7772499e+007 
Iter:  9 Ap: 2.59e-001 Pobj: -1.7872294e-001 Ad: 1.16e-001 Dobj: -2.6204276e+007 
Iter: 10 Ap: 3.30e-001 Pobj: -1.4676353e-001 Ad: 3.77e-001 Dobj: -2.3152780e+007 
Iter: 11 Ap: 2.68e-001 Pobj: -8.8779044e-002 Ad: 2.48e-001 Dobj: -2.6689109e+007 
Iter: 12 Ap: 1.96e-001 Pobj: -7.1291821e-002 Ad: 1.73e-001 Dobj: -3.3215252e+007 
Iter: 13 Ap: 1.72e-001 Pobj: -8.8379535e-002 Ad: 1.35e-001 Dobj: -3.6279903e+007 
Iter: 14 Ap: 1.64e-001 Pobj: -6.8449424e-002 Ad: 1.52e-001 Dobj: -4.7157744e+007 
Iter: 15 Ap: 2.37e-001 Pobj: -3.1114869e-002 Ad: 1.89e-001 Dobj: -4.6955571e+007 
Iter: 16 Ap: 1.71e-001 Pobj: -2.2216149e-002 Ad: 1.02e-001 Dobj: -4.9721674e+007 
Iter: 17 Ap: 1.14e-001 Pobj: -2.8719540e-002 Ad: 1.35e-001 Dobj: -6.3235837e+007 
Iter: 18 Ap: 1.73e-001 Pobj: -2.2907474e-002 Ad: 8.54e-002 Dobj: -7.0021527e+007 
Iter: 19 Ap: 1.56e-001 Pobj: -6.3163150e-002 Ad: 1.37e-001 Dobj: -9.1592672e+007 
Iter: 20 Ap: 1.19e-001 Pobj: -8.4291142e-002 Ad: 9.67e-002 Dobj: -1.1694313e+008 
Iter: 21 Ap: 6.65e-002 Pobj: -9.0550316e-002 Ad: 7.97e-002 Dobj: -1.4573091e+008 
Iter: 22 Ap: 1.07e-001 Pobj: -1.0265732e-001 Ad: 7.92e-002 Dobj: -1.4335241e+008 
Iter: 23 Ap: 1.40e-001 Pobj: -1.2893908e-001 Ad: 9.89e-002 Dobj: -1.6888717e+008 
Iter: 24 Ap: 1.04e-001 Pobj: -1.2968038e-001 Ad: 8.45e-002 Dobj: -2.0106444e+008 
Iter: 25 Ap: 1.19e-001 Pobj: -1.3396755e-001 Ad: 1.03e-001 Dobj: -2.6291206e+008 
Iter: 26 Ap: 1.61e-001 Pobj: -1.2081859e-001 Ad: 1.29e-001 Dobj: -3.1272945e+008 
Iter: 27 Ap: 1.47e-001 Pobj: -1.6442645e-001 Ad: 6.79e-002 Dobj: -3.4553305e+008 
Iter: 28 Ap: 1.90e-001 Pobj: -1.7087796e-001 Ad: 1.83e-001 Dobj: -4.1812515e+008 
Iter: 29 Ap: 1.46e-001 Pobj: -1.6625690e-001 Ad: 1.19e-001 Dobj: -4.9763111e+008 
Iter: 30 Ap: 1.79e-001 Pobj: -1.5523266e-001 Ad: 1.11e-001 Dobj: -5.3710511e+008 
Iter: 31 Ap: 1.05e-001 Pobj: -1.3872120e-001 Ad: 5.70e-002 Dobj: -5.8421301e+008 
Iter: 32 Ap: 2.46e-001 Pobj: -1.2280831e-001 Ad: 2.46e-001 Dobj: -6.4857491e+008 
Iter: 33 Ap: 2.04e-001 Pobj: -1.1117251e-001 Ad: 2.55e-001 Dobj: -7.0580803e+008 
Iter: 34 Ap: 2.46e-001 Pobj: -1.1812814e-001 Ad: 3.40e-001 Dobj: -7.0138268e+008 
Iter: 35 Ap: 1.70e-001 Pobj: -1.0929292e-001 Ad: 1.95e-001 Dobj: -6.7625417e+008 
Iter: 36 Ap: 8.32e-002 Pobj: -9.1785166e-002 Ad: 9.16e-002 Dobj: -6.5604306e+008 
Iter: 37 Ap: 1.28e-001 Pobj: -8.1225888e-002 Ad: 1.50e-001 Dobj: -6.2504958e+008 
Iter: 38 Ap: 2.26e-001 Pobj: -6.8912568e-002 Ad: 2.25e-001 Dobj: -5.9836368e+008 
Iter: 39 Ap: 3.10e-001 Pobj: -5.8786865e-002 Ad: 2.12e-001 Dobj: -5.7828546e+008 
Iter: 40 Ap: 9.61e-002 Pobj: -1.2469784e-002 Ad: 3.67e-002 Dobj: -5.7263015e+008 
Iter: 41 Ap: 4.23e-001 Pobj:  1.0058632e-001 Ad: 1.98e-001 Dobj: -5.4950178e+008 
Iter: 42 Ap: 7.45e-002 Pobj:  2.4967065e-001 Ad: 8.00e-002 Dobj: -5.3450070e+008 
Iter: 43 Ap: 5.48e-002 Pobj:  3.8282919e-001 Ad: 8.20e-002 Dobj: -5.1029128e+008 
Iter: 44 Ap: 1.74e-001 Pobj:  5.9074755e-001 Ad: 1.16e-001 Dobj: -4.9590747e+008 
Iter: 45 Ap: 9.47e-002 Pobj:  1.0574616e+000 Ad: 1.32e-001 Dobj: -4.7476221e+008 
Iter: 46 Ap: 1.36e-001 Pobj:  1.3202404e+000 Ad: 1.28e-001 Dobj: -4.6041614e+008 
Iter: 47 Ap: 2.89e-001 Pobj:  6.8435726e+000 Ad: 8.00e-002 Dobj: -4.5444093e+008 
Iter: 48 Ap: 1.67e-001 Pobj:  8.7056222e+000 Ad: 1.49e-001 Dobj: -4.3147855e+008 
Iter: 49 Ap: 3.75e-002 Pobj:  1.0702630e+001 Ad: 7.64e-002 Dobj: -4.1525341e+008 
Iter: 50 Ap: 1.21e-001 Pobj:  2.4260114e+001 Ad: 1.33e-001 Dobj: -4.0489770e+008 
Iter: 51 Ap: 4.95e-002 Pobj:  4.5650516e+001 Ad: 2.62e-002 Dobj: -4.0126904e+008 
Iter: 52 Ap: 1.89e-002 Pobj:  5.5841017e+001 Ad: 6.13e-002 Dobj: -3.9325771e+008 
Iter: 53 Ap: 2.98e-003 Pobj:  5.7611261e+001 Ad: 5.44e-002 Dobj: -3.8728660e+008 
Iter: 54 Ap: 1.28e-003 Pobj:  5.8058111e+001 Ad: 7.24e-002 Dobj: -3.7893624e+008 
Iter: 55 Ap: 6.69e-003 Pobj:  6.2474520e+001 Ad: 2.34e-002 Dobj: -3.7531119e+008 
Iter: 56 Ap: 2.39e-003 Pobj:  6.3323565e+001 Ad: 3.04e-002 Dobj: -3.7239800e+008 
Iter: 57 Ap: 7.90e-003 Pobj:  7.0040754e+001 Ad: 5.02e-002 Dobj: -3.6727162e+008 
Iter: 58 Ap: 4.98e-004 Pobj:  7.0251734e+001 Ad: 6.70e-002 Dobj: -3.6011516e+008 
Iter: 59 Ap: 7.13e-003 Pobj:  7.4779127e+001 Ad: 2.68e-002 Dobj: -3.5509648e+008 
Iter: 60 Ap: 2.02e-003 Pobj:  7.6522387e+001 Ad: 6.77e-002 Dobj: -3.5012660e+008 
Iter: 61 Ap: 2.38e-003 Pobj:  7.9244340e+001 Ad: 3.61e-002 Dobj: -3.4788780e+008 
Iter: 62 Ap: 5.54e-003 Pobj:  8.6349968e+001 Ad: 2.54e-002 Dobj: -3.4479584e+008 
Iter: 63 Ap: 1.24e-003 Pobj:  8.8030114e+001 Ad: 5.35e-002 Dobj: -3.4137316e+008 
Iter: 64 Ap: 2.84e-003 Pobj:  9.1882687e+001 Ad: 3.62e-002 Dobj: -3.3848881e+008 
Lack of progress.  Giving up!
Failure: return code is 7 
Primal objective value: 1.0702630e+001
Dual objective value: -4.1525341e+008
Relative primal infeasibility: 2.76e-001
Relative dual infeasibility: 5.31e+001
Real Relative Gap: -1.00e+000
XZ Relative Gap: 4.25e+001
DIMACS error measures: 3.94e-001 0.00e+000 9.17e+001 0.00e+000 -1.00e+000 4.25e+001
```

那么无论如何，先把剩下的测试做完吧。

运行诊断代码
```julia
using JLD2
working_path = "D:\\Projects\\numerical-boostrap\\oscillator-simple-prototype\\"
@load working_path * "xpopstr_expected_ode-dx-0.02-l-10.jld2" xpopstr_expected_ode 

xpopstr_expected_ope_value = zeros(2xpopspace_dim)

for xpopstr_idx in 1 : xpopspace_dim
    x_power = index_to_xpower(xpopstr_idx)
    p_power = index_to_ppower(xpopstr_idx)
    expected_value = xpopstr_expected_ode[x_power, p_power]
    xpopstr_expected_ope_value[2xpopstr_idx - 1] = real(expected_value)
    xpopstr_expected_ope_value[2xpopstr_idx] = imag(expected_value)
end

xpopstr_expected_real_imag_parts_ope_value(i, real_or_imag) = begin
    if real_or_imag == :real
        return xpopstr_expected_ope_value[2i - 1]
    end
    if real_or_imag == :imag
        return xpopstr_expected_ope_value[2i]
    end
end

variable_list_real_ope_value = 
    [xpopstr_expected_real_imag_parts_ope_value(i, :real) * I22 for i in 1 : xpopspace_dim]
variable_list_imag_ope_value = 
    [xpopstr_expected_real_imag_parts_ope_value(i, :imag) * Im22 for i in 1 : xpopspace_dim]
xpopstr_basis_real_ope_value = OffsetArray([I22, variable_list_real_ope_value...], xpopspace_index_range)
xpopstr_basis_imag_ope_value = OffsetArray([O22, variable_list_imag_ope_value...], xpopspace_index_range)

function complex_to_mat_ope_value(coefficients)
    real_part = transpose(real(coefficients))
    imag_part = transpose(imag(coefficients))
    real_part_mat_version = map(x -> x * I22, real_part)
    imag_part_mat_version = map(x -> x * Im22, imag_part)
    (real_part_mat_version + imag_part_mat_version) * (xpopstr_basis_real_ope_value + xpopstr_basis_imag_ope_value)
end

for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = real(cons)
    cons_imag = imag(cons)
    lhs = complex_to_mat_ope_value(cons)
    if cons_real == cons_imag == zero_xpopstr
        continue
    end
    if cons_real == zero_xpopstr
        println("$x_power $p_power   $(lhs[1, 2])")
    elseif cons_imag == zero_xpopstr
        println("$x_power $p_power   $(lhs[1, 1])")
    else
        println("$x_power $p_power   $(lhs[1, 2])")
        println("$x_power $p_power   $(lhs[1, 1])")
    end
end
```
发现
```
0 1   -4.902451214784239e-14
0 2   0.0
0 2   -0.0006943306158411211
0 3   -0.32820777442942706
0 3   0.0
1 0   0.0
1 1   0.0006444451709375354
1 2   0.0
1 2   -0.17131613487814604
1 3   0.020025042775518642
1 3   0.0
```

我们下面讨论
```
0 3   -0.32820777442942706
```
这一行。计算
```julia
xpopstr_stringify(comm_with_ham(xpopstr_xp_power(0, 3)))
```
得到
```
- 6.0 p + 24.0im x - 6.0im x p^2 - 36.0 x^2 p - 12.0im x^3 p^2
```
这个点看起来是纯粹的数值误差。我们有
```
julia> xpopstr_expected_ode[0, 1], xpopstr_expected_ode[1, 0], xpopstr_expected_ode[1, 2], xpopstr_expected_ode[2, 1], xpopstr_expected_ode[3, 2]
(0.0 + 0.008262813010199016im, -9.672537983335926e-15 + 0.0im, -0.024780705457444523 + 0.0im, 0.0 - 0.0018168737481579495im, -0.013641080401041665 + 0.0im)
```
因此应该是计算误差积累的结果。

不过我倒是意外发现了一个bug，就是
```julia
for x_power in 0 : L_max - 4, p_power in 0 : L_max - 2
```
似乎应该替换成
```julia
for x_power in 0 : 2L_max - 4, p_power in 0 : 2L_max - 2
```

总结一下目前的进展吧。$\lang O, H \rang = 0$的bug可能没有多少了，正定性的bug也应该没有多少了，为了判断是否真是如此，拿Mathematica算benchmark将是当务之急；如果事后发现Mathematica给出的解确实是feasible的，那么应该说主要问题在优化的技术细节上面了。

此外，改用`COSMO`以后，目标函数变成了-20。可能是SDP约束不到位？

运行诊断代码
```julia
using JLD2
working_path = "D:\\Projects\\numerical-boostrap\\oscillator-simple-prototype\\"
@load working_path * "xpopstr_expected_ode-dx-0.02-l-10.jld2" xpopstr_expected_ode 

xpopstr_expected_ope_value = zeros(2xpopspace_dim)

for xpopstr_idx in 1 : xpopspace_dim
    x_power = index_to_xpower(xpopstr_idx)
    p_power = index_to_ppower(xpopstr_idx)
    expected_value = xpopstr_expected_ode[x_power, p_power]
    xpopstr_expected_ope_value[2xpopstr_idx - 1] = real(expected_value)
    xpopstr_expected_ope_value[2xpopstr_idx] = imag(expected_value)
end

xpopstr_expected_real_imag_parts_ope_value(i, real_or_imag) = begin
    if real_or_imag == :real
        return xpopstr_expected_ope_value[2i - 1]
    end
    if real_or_imag == :imag
        return xpopstr_expected_ope_value[2i]
    end
end

variable_list_real_ope_value = 
    [xpopstr_expected_real_imag_parts_ope_value(i, :real) * I22 for i in 1 : xpopspace_dim]
variable_list_imag_ope_value = 
    [xpopstr_expected_real_imag_parts_ope_value(i, :imag) * Im22 for i in 1 : xpopspace_dim]
xpopstr_basis_real_ope_value = OffsetArray([I22, variable_list_real_ope_value...], xpopspace_index_range)
xpopstr_basis_imag_ope_value = OffsetArray([O22, variable_list_imag_ope_value...], xpopspace_index_range)

function complex_to_mat_ope_value(coefficients)
    real_part = transpose(real(coefficients))
    imag_part = transpose(imag(coefficients))
    real_part_mat_version = map(x -> x * I22, real_part)
    imag_part_mat_version = map(x -> x * Im22, imag_part)
    (real_part_mat_version + imag_part_mat_version) * (xpopstr_basis_real_ope_value + xpopstr_basis_imag_ope_value)
end

for x_power in 0 : 2L_max - 4, p_power in 0 : 2L_max - 2
    op = xpopstr_xp_power(x_power, p_power)
    cons = comm_with_ham(op)
    cons_real = real(cons)
    cons_imag = imag(cons)
    lhs = complex_to_mat_ope_value(cons)
    if cons_real == cons_imag == zero_xpopstr
        continue
    end
    if cons_real == zero_xpopstr
        println("$x_power $p_power   $(lhs[1, 2])")
    elseif cons_imag == zero_xpopstr
        println("$x_power $p_power   $(lhs[1, 1])")
    else
        println("$x_power $p_power   $(lhs[1, 2])")
        println("$x_power $p_power   $(lhs[1, 1])")
    end
end

##

M_ode_value = zeros(2 * (L_max + 1)^2, 2 * (L_max + 1)^2)

for i in 1 : (L_max + 1)^2
    for j in i : (L_max + 1)^2
        op1_idx = M_index_to_xpopstr_index[i]
        op2_idx = M_index_to_xpopstr_index[j]
        op1_idx_xpower = index_to_xpower(op1_idx)
        op1_idx_ppower = index_to_ppower(op1_idx)
        op2_idx_xpower = index_to_xpower(op2_idx)
        op2_idx_ppower = index_to_ppower(op2_idx)
        op_ij = xpopstr_normal_ord(op1_idx_xpower, op1_idx_ppower, op2_idx_xpower, op2_idx_ppower)

        M_ode_value[2i - 1 : 2i, 2j - 1 : 2j] = complex_to_mat_ope_value(op_ij)
    end
end

plot(eigen(M_ode_value).values, legend=false)
```
得到的结果见`m-eigen-value-2022-3-20.png`。无非两种情况：
- 要么解方程得到的关联函数值非常不精确（目前的配置：`Δx = 0.02; L = 10`）
- 要么约束就写错了

## 2022.3.21

在`calculate-all-correlation-functions-in-xp-oscillator-for-benchmark-3.nb`中有存放$\lang x^n p^m \rang$的数据。计算策略在这个笔记本里面写得很清楚了。
计算时，首先使用$x^n, n \leq 8$的nonlinear SDP做计算。所得非零部分列举如下：
```Mathematica
{xpOpString[] -> 1, xpOpString[p, p] -> 0.808702 + 0. I, 
 xpOpString[x, p] -> 0. + 0.5 I, xpOpString[x, x] -> 0.301138, 
 xpOpString[p, p, p, p] -> 1.88474 + 0. I, 
 xpOpString[x, p, p, p] -> 0. + 1.21305 I, 
 xpOpString[x, x, p, p] -> -0.186501 + 0. I, 
 xpOpString[x, x, x, p] -> 0. + 0.451707 I, 
 xpOpString[x, x, x, x] -> 0.253782, 
 xpOpString[p, p, p, p, p, p] -> 7.01002 + 0. I, 
 xpOpString[x, p, p, p, p, p] -> 0. + 4.71184 I, 
 xpOpString[x, x, p, p, p, p] -> -1.45084 + 0. I, 
 xpOpString[x, x, x, p, p, p] -> 0. + 0.660746 I, 
 xpOpString[x, x, x, x, p, p] -> -0.595472 + 0. I, 
 xpOpString[x, x, x, x, x, p] -> 0. + 0.634454 I, 
 xpOpString[x, x, x, x, x, x] -> 0.343358, 
 xpOpString[p, p, p, p, p, p, p, p] -> 36.8242 + 0. I, 
 xpOpString[x, p, p, p, p, p, p, p] -> 0. + 24.5351 I, 
 xpOpString[x, x, p, p, p, p, p, p] -> -8.90225 + 0. I, 
 xpOpString[x, x, x, p, p, p, p, p] -> 0. + 1.24922 I, 
 xpOpString[x, x, x, x, p, p, p, p] -> -2.98912 + 0. I, 
 xpOpString[x, x, x, x, x, p, p, p] -> 0. + 0.051031 I, 
 xpOpString[x, x, x, x, x, x, p, p] -> -1.44281 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, p] -> 0. + 1.20175 I, 
 xpOpString[x, x, x, x, x, x, x, x] -> 0.598177, 
 xpOpString[p, p, p, p, p, p, p, p, p, p] -> 280.41 + 0. I, 
 xpOpString[x, p, p, p, p, p, p, p, p, p] -> 0. + 165.709 I, 
 xpOpString[x, x, p, p, p, p, p, p, p, p] -> -56.5456 + 0. I, 
 xpOpString[x, x, x, p, p, p, p, p, p, p] -> 0. + 5.47501 I, 
 xpOpString[x, x, x, x, p, p, p, p, p, p] -> -18.363 + 0. I, 
 xpOpString[x, x, x, x, x, p, p, p, p, p] -> 0. - 5.33909 I, 
 xpOpString[x, x, x, x, x, x, p, p, p, p] -> -5.16492 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, p, p, p] -> 0. - 1.82597 I, 
 xpOpString[x, x, x, x, x, x, x, x, p, p] -> -3.90405 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p] -> 0. + 2.6918 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x] -> 1.31284, 
 xpOpString[x, x, p, p, p, p, p, p, p, p, p, p] -> -378.595 + 0. I, 
 xpOpString[x, x, x, p, p, p, p, p, p, p, p, p] -> 0. + 119.898 I, 
 xpOpString[x, x, x, x, p, p, p, p, p, p, p, p] -> -173.651 + 0. I, 
 xpOpString[x, x, x, x, x, p, p, p, p, p, p, p] -> 0. - 64.0791 I, 
 xpOpString[x, x, x, x, x, x, p, p, p, p, p, p] -> -21.351 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, p, p, p, p, p] -> 0. - 23.5749 I, 
 xpOpString[x, x, x, x, x, x, x, x, p, p, p, p] -> -7.952 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p, p, p] -> 0. - 9.44159 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, p, p] -> -11.375 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, p] -> 0. + 7.22064 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x] -> 3.40689, 
 xpOpString[x, x, x, x, p, p, p, p, p, p, p, p, p, p] -> -2545.91 + 
   0. I, xpOpString[x, x, x, x, x, p, p, p, p, p, p, p, p, p] -> 
  0. - 875.37 I, 
 xpOpString[x, x, x, x, x, x, p, p, p, p, p, p, p, p] -> -114.049 + 
   0. I, xpOpString[x, x, x, x, x, x, x, p, p, p, p, p, p, p] -> 
  0. - 240.416 I, 
 xpOpString[x, x, x, x, x, x, x, x, p, p, p, p, p, p] -> 
  19.3344 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p, p, p, p, p] -> 
  0. - 78.2718 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, p, p, p, p] -> -3.34435 + 
   0. I, xpOpString[x, x, x, x, x, x, x, x, x, x, x, p, p, p] -> 
  0. - 39.639 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, p, p] -> -37.5785 + 
   0. I, xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, p] -> 
  0. + 22.1448 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x] -> 9.75835, 
 xpOpString[x, x, x, x, x, x, p, p, p, p, p, p, p, p, p, 
   p] -> -718.932 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, p, p, p, p, p, p, p, p, p] -> 
  0. - 3105.63 I, 
 xpOpString[x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p] -> 
  785.441 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p] -> 
  0. - 598.6 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p] -> 
  283.659 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p] -> 
  0. - 236.616 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p] -> 
  61.1296 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p] -> 
  0. - 169.57 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, 
   p] -> -137.121 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p] -> 
  0. + 73.1876 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x] -> 
  32.4658, xpOpString[x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p, 
   p, p] -> 17735.1 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p, p] -> 
  0. - 5504.13 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p] -> 
  5072.89 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p] -> 
  0. - 911.61 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p] -> 
  1608.26 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p] -> 
  0. - 620.892 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p] -> 
  500.407 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p] -> 
  0. - 760.017 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, 
   p] -> -525.418 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p] -> 
  0. + 275.96 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x] -> 
  117.962, xpOpString[x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, 
   p, p, p, p] -> 93627.8 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p, 
   p] -> 0. + 7321.42 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, 
   p] -> 21566.4 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, 
   p] -> 0. + 2427.49 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, 
   p] -> 7417.75 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, 
   p] -> 0. - 1159.55 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, 
   p] -> 3212.02 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, 
   p] -> 0. - 3444.63 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, 
   p] -> -2253.29 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, 
   p] -> 0. + 1120.64 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, 
   x] -> 451.727}
```
由此获得的$M$矩阵好像并不正定……还有一个问题是，似乎$\lang x^2 p^2 \rang$是负的。哦不过这个好像确实是对的，之前解方程的时候算出来过。
特征值是
```
{-72560.2 + 5.22106*10^-19 I, -18647.5 - 3.62637*10^-14 I, -3472.86 + 
  3.88545*10^-13 I, -2151.5 + 2.06363*10^-18 I, -184.817 + 
  3.37672*10^-17 I, 92.3838 - 77.3257 I, 
 92.3838 + 77.3257 I, -41.2894 + 5.70591*10^-16 I, -19.772 - 
  1.35127*10^-15 I, -6.29202 - 2.62852*10^-14 I, 2.75494 - 3.62274 I, 
 2.75494 + 3.62274 I, -3.72358 + 6.43788*10^-13 I, 
 1.95103 + 1.27914 I, 1.95103 - 1.27914 I, 0.343084 - 0.295659 I, 
 0.343084 + 0.295659 I, -0.302859 + 0.0292939 I, -0.302859 - 
  0.0292939 I, -0.0941271 + 1.48248*10^-13 I, -0.00625381 + 
  0.0170569 I, -0.00625381 - 0.0170569 I, -0.0153286 - 
  1.13102*10^-13 I, 
 0.0100093 - 2.65962*10^-13 I, -0.00182255 - 
  0.00680812 I, -0.00182255 + 0.00680812 I, -0.00622958 + 
  1.45031*10^-13 I, 
 0.00338418 - 3.70776*10^-14 I, -1.53504*10^-12 + 9.65472*10^-13 I, 
 1.25553*10^-12 + 1.03573*10^-12 I, 
 6.45246*10^-13 + 6.56672*10^-13 I, -3.19067*10^-13 - 
  7.00729*10^-13 I, -3.16211*10^-13 - 3.92248*10^-13 I, 
 3.86535*10^-13 - 9.19704*10^-15 I, 2.32805*10^-13 - 2.83706*10^-13 I,
  4.47278*10^-14 + 1.03133*10^-13 I}
```
另一方面，我们有
```Mathematica
Product[Eigenvalues[(M /. nonzeroExpected /. 
     xpOpString[seq__] /; OddQ[Length[{seq}]] -> 0)][[i]], {i, 1, 
  Length[M]}]
```
输出
```
9.29718*10^-89 + 6.60311*10^-90 I
```

实际上，如果我们不包含那么高次的算符，就有
```Mathematica
Eigenvalues[(Table[
     xpNormalOrderedOp[i + j, 0], {i, 0, 10}, {j, 0, 10}] /. 
    nonzeroExpected /. xpOpString[seq__] /; OddQ[Length[{seq}]] -> 0)]
```
输出
```
485.275, 127.919, 3.1433, 1.58921, 0.955657, 0.227242, 0.0831391, \
-0.0774592, -0.0658394, 0.0608265, 0.0198087
```

最后，或许可以考虑一下检查feasibility。能够通过增大容差来获得勉强令人满意的结果？

以下是$x^n, n\leq 11$的结果：
```
{xpOpString[] -> 1, xpOpString[p, p] -> 0.826205 + 0. I, 
 xpOpString[x, p] -> 0. + 0.5 I, xpOpString[x, x] -> 0.305782, 
 xpOpString[p, p, p, p] -> 1.93335 + 0. I, 
 xpOpString[x, p, p, p] -> 0. + 1.23931 I, 
 xpOpString[x, x, p, p] -> -0.181759 + 0. I, 
 xpOpString[x, x, x, p] -> 0. + 0.458672 I, 
 xpOpString[x, x, x, x] -> 0.260212, 
 xpOpString[p, p, p, p, p, p] -> 7.2212 + 0. I, 
 xpOpString[x, p, p, p, p, p] -> 0. + 4.83338 I, 
 xpOpString[x, x, p, p, p, p] -> -1.47756 + 0. I, 
 xpOpString[x, x, x, p, p, p] -> 0. + 0.682085 I, 
 xpOpString[x, x, x, x, p, p] -> -0.601349 + 0. I, 
 xpOpString[x, x, x, x, x, p] -> 0. + 0.65053 I, 
 xpOpString[x, x, x, x, x, x] -> 0.347256, 
 xpOpString[p, p, p, p, p, p, p, p] -> 38.0874 + 0. I, 
 xpOpString[x, p, p, p, p, p, p, p] -> 0. + 25.2742 I, 
 xpOpString[x, x, p, p, p, p, p, p] -> -9.10923 + 0. I, 
 xpOpString[x, x, x, p, p, p, p, p] -> 0. + 1.31136 I, 
 xpOpString[x, x, x, x, p, p, p, p] -> -3.05202 + 0. I, 
 xpOpString[x, x, x, x, x, p, p, p] -> 0. + 0.0766033 I, 
 xpOpString[x, x, x, x, x, x, p, p] -> -1.47895 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, p] -> 0. + 1.21539 I, 
 xpOpString[x, x, x, x, x, x, x, x] -> 0.61636, 
 xpOpString[p, p, p, p, p, p, p, p, p, p] -> 290.984 + 0. I, 
 xpOpString[x, p, p, p, p, p, p, p, p, p] -> 0. + 171.393 I, 
 xpOpString[x, x, p, p, p, p, p, p, p, p] -> -58.4308 + 0. I, 
 xpOpString[x, x, x, p, p, p, p, p, p, p] -> 0. + 5.85395 I, 
 xpOpString[x, x, x, x, p, p, p, p, p, p] -> -18.9248 + 0. I, 
 xpOpString[x, x, x, x, x, p, p, p, p, p] -> 0. - 5.41413 I, 
 xpOpString[x, x, x, x, x, x, p, p, p, p] -> -5.36261 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, p, p, p] -> 0. - 1.86789 I, 
 xpOpString[x, x, x, x, x, x, x, x, p, p] -> -3.94401 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p] -> 0. + 2.77362 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x] -> 1.34604, 
 xpOpString[x, x, p, p, p, p, p, p, p, p, p, p] -> -393.285 + 0. I, 
 xpOpString[x, x, x, p, p, p, p, p, p, p, p, p] -> 0. + 121.055 I, 
 xpOpString[x, x, x, x, p, p, p, p, p, p, p, p] -> -179.528 + 0. I, 
 xpOpString[x, x, x, x, x, p, p, p, p, p, p, p] -> 0. - 65.8857 I, 
 xpOpString[x, x, x, x, x, x, p, p, p, p, p, p] -> -21.9806 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, p, p, p, p, p] -> 0. - 24.2693 I, 
 xpOpString[x, x, x, x, x, x, x, x, p, p, p, p] -> -8.18325 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p, p, p] -> 0. - 9.48987 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, p, p] -> -11.7121 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, p] -> 0. + 7.40324 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x] -> 3.45606, 
 xpOpString[x, x, x, x, p, p, p, p, p, p, p, p, p, p] -> -2602.5 + 
   0. I, xpOpString[x, x, x, x, x, p, p, p, p, p, p, p, p, p] -> 
  0. - 900.886 I, 
 xpOpString[x, x, x, x, x, x, p, p, p, p, p, p, p, p] -> -117.932 + 
   0. I, xpOpString[x, x, x, x, x, x, x, p, p, p, p, p, p, p] -> 
  0. - 245.959 I, 
 xpOpString[x, x, x, x, x, x, x, x, p, p, p, p, p, p] -> 
  19.0191 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p, p, p, p, p] -> 
  0. - 80.4034 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, p, p, p, p] -> -3.89309 + 
   0. I, xpOpString[x, x, x, x, x, x, x, x, x, x, x, p, p, p] -> 
  0. - 40.7005 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, p, p] -> -38.5306 + 
   0. I, xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, p] -> 
  0. + 22.4644 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x] -> 10.13, 
 xpOpString[x, x, x, x, x, x, p, p, p, p, p, p, p, p, p, 
   p] -> -780.606 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, p, p, p, p, p, p, p, p, p] -> 
  0. - 3199.36 I, 
 xpOpString[x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p] -> 
  800.881 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p] -> 
  0. - 623.083 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p] -> 
  291.196 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p] -> 
  0. - 242.547 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p] -> 
  61.4837 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p] -> 
  0. - 173.895 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, 
   p] -> -139.045 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p] -> 
  0. + 75.975 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x] -> 
  33.2121, xpOpString[x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p, 
   p, p] -> 18280.6 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p, p] -> 
  0. - 5645.39 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p] -> 
  5240.32 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p] -> 
  0. - 938.236 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p] -> 
  1638.83 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p] -> 
  0. - 651.351 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p] -> 
  513.476 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p] -> 
  0. - 769.754 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, 
   p] -> -545.267 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p] -> 
  0. + 282.303 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x] -> 
  119.936, xpOpString[x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, 
   p, p, p, p] -> 95736.1 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, p, 
   p] -> 0. + 7266.96 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, p, 
   p] -> 22179.5 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, p, 
   p] -> 0. + 2422.99 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, p, 
   p] -> 7699.63 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, p, 
   p] -> 0. - 1186.15 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, p, 
   p] -> 3240.54 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, p, 
   p] -> 0. - 3571.69 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, p, 
   p] -> -2305.31 + 0. I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, 
   p] -> 0. + 1139.39 I, 
 xpOpString[x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, x, 
   x] -> 471.569}
```
正定性的情况有所好转：
```Mathematica
Eigenvalues[(Table[
     xpNormalOrderedOp[i + j, 0], {i, 0, 10}, {j, 0, 10}] /. 
    nonzeroExpected /. xpOpString[seq__] /; OddQ[Length[{seq}]] -> 0)]
```
给出
```
{504.899, 130.211, 4.0255, 1.59372, 1.01345, 0.236806, 0.16176, \
0.0226981, 0.0132164, 0.000721841, 0.000523045}
```
此时
```Mathematica
Product[Eigenvalues[(M /. nonzeroExpected /. 
     xpOpString[seq__] /; OddQ[Length[{seq}]] -> 0)][[i]], {i, 1, 
  Length[M]}]
```
给出
```
-1.28092*10^-123 + 1.06183*10^-123 I
```

此时的$M$矩阵保存在`nonlinear-SDP-x-power-10-standard-m-value.csv`中。

现在`xpopstr_expected`和`M`都有了，那么可以检查`jump-oscillator-2-benchmark-with-ode.jl`是不是正确了。

在`obtain-jump-point-xpopstr-expected-nonlinear-sdp-x-power-11.nb`中自动生成用于JuMP feasibility测试的数据点：
```julia
xpopstr_expected[1] => 0, xpopstr_expected[2] => 0, 
xpopstr_expected[3] => 0.826205, xpopstr_expected[4] => 0., 
xpopstr_expected[5] => 0, xpopstr_expected[6] => 0, 
xpopstr_expected[7] => 1.93335, xpopstr_expected[8] => 0., 
xpopstr_expected[9] => 0, xpopstr_expected[10] => 0, 
xpopstr_expected[11] => 7.2212, xpopstr_expected[12] => 0., 
xpopstr_expected[13] => 0, xpopstr_expected[14] => 0, 
xpopstr_expected[15] => 38.0874, xpopstr_expected[16] => 0., 
xpopstr_expected[17] => 0, xpopstr_expected[18] => 0, 
xpopstr_expected[19] => 290.984, xpopstr_expected[20] => 0., 
xpopstr_expected[21] => 0, xpopstr_expected[22] => 0, 
xpopstr_expected[23] => 0., xpopstr_expected[24] => 0.5, 
xpopstr_expected[25] => 0, xpopstr_expected[26] => 0, 
xpopstr_expected[27] => 0., xpopstr_expected[28] => 1.23931, 
xpopstr_expected[29] => 0, xpopstr_expected[30] => 0, 
xpopstr_expected[31] => 0., xpopstr_expected[32] => 4.83338, 
xpopstr_expected[33] => 0, xpopstr_expected[34] => 0, 
xpopstr_expected[35] => 0., xpopstr_expected[36] => 25.2742, 
xpopstr_expected[37] => 0, xpopstr_expected[38] => 0, 
xpopstr_expected[39] => 0., xpopstr_expected[40] => 171.393, 
xpopstr_expected[41] => 0, xpopstr_expected[42] => 0, 
xpopstr_expected[43] => 0.305782, xpopstr_expected[44] => 0, 
xpopstr_expected[45] => 0, xpopstr_expected[46] => 0, 
xpopstr_expected[47] => -0.181759, xpopstr_expected[48] => 0., 
xpopstr_expected[49] => 0, xpopstr_expected[50] => 0, 
xpopstr_expected[51] => -1.47756, xpopstr_expected[52] => 0., 
xpopstr_expected[53] => 0, xpopstr_expected[54] => 0, 
xpopstr_expected[55] => -9.10923, xpopstr_expected[56] => 0., 
xpopstr_expected[57] => 0, xpopstr_expected[58] => 0, 
xpopstr_expected[59] => -58.4308, xpopstr_expected[60] => 0., 
xpopstr_expected[61] => 0, xpopstr_expected[62] => 0, 
xpopstr_expected[63] => -393.285, xpopstr_expected[64] => 0., 
xpopstr_expected[65] => 0, xpopstr_expected[66] => 0, 
xpopstr_expected[67] => 0., xpopstr_expected[68] => 0.458672, 
xpopstr_expected[69] => 0, xpopstr_expected[70] => 0, 
xpopstr_expected[71] => 0., xpopstr_expected[72] => 0.682085, 
xpopstr_expected[73] => 0, xpopstr_expected[74] => 0, 
xpopstr_expected[75] => 0., xpopstr_expected[76] => 1.31136, 
xpopstr_expected[77] => 0, xpopstr_expected[78] => 0, 
xpopstr_expected[79] => 0., xpopstr_expected[80] => 5.85395, 
xpopstr_expected[81] => 0, xpopstr_expected[82] => 0, 
xpopstr_expected[83] => 0., xpopstr_expected[84] => 121.055, 
xpopstr_expected[85] => 0, xpopstr_expected[86] => 0, 
xpopstr_expected[87] => 0.260212, xpopstr_expected[88] => 0, 
xpopstr_expected[89] => 0, xpopstr_expected[90] => 0, 
xpopstr_expected[91] => -0.601349, xpopstr_expected[92] => 0., 
xpopstr_expected[93] => 0, xpopstr_expected[94] => 0, 
xpopstr_expected[95] => -3.05202, xpopstr_expected[96] => 0., 
xpopstr_expected[97] => 0, xpopstr_expected[98] => 0, 
xpopstr_expected[99] => -18.9248, xpopstr_expected[100] => 0., 
xpopstr_expected[101] => 0, xpopstr_expected[102] => 0, 
xpopstr_expected[103] => -179.528, xpopstr_expected[104] => 0., 
xpopstr_expected[105] => 0, xpopstr_expected[106] => 0, 
xpopstr_expected[107] => -2602.5, xpopstr_expected[108] => 0., 
xpopstr_expected[109] => 0, xpopstr_expected[110] => 0, 
xpopstr_expected[111] => 0., xpopstr_expected[112] => 0.65053, 
xpopstr_expected[113] => 0, xpopstr_expected[114] => 0, 
xpopstr_expected[115] => 0., xpopstr_expected[116] => 0.0766033, 
xpopstr_expected[117] => 0, xpopstr_expected[118] => 0, 
xpopstr_expected[119] => 0., xpopstr_expected[120] => -5.41413, 
xpopstr_expected[121] => 0, xpopstr_expected[122] => 0, 
xpopstr_expected[123] => 0., xpopstr_expected[124] => -65.8857, 
xpopstr_expected[125] => 0, xpopstr_expected[126] => 0, 
xpopstr_expected[127] => 0., xpopstr_expected[128] => -900.886, 
xpopstr_expected[129] => 0, xpopstr_expected[130] => 0, 
xpopstr_expected[131] => 0.347256, xpopstr_expected[132] => 0, 
xpopstr_expected[133] => 0, xpopstr_expected[134] => 0, 
xpopstr_expected[135] => -1.47895, xpopstr_expected[136] => 0., 
xpopstr_expected[137] => 0, xpopstr_expected[138] => 0, 
xpopstr_expected[139] => -5.36261, xpopstr_expected[140] => 0., 
xpopstr_expected[141] => 0, xpopstr_expected[142] => 0, 
xpopstr_expected[143] => -21.9806, xpopstr_expected[144] => 0., 
xpopstr_expected[145] => 0, xpopstr_expected[146] => 0, 
xpopstr_expected[147] => -117.932, xpopstr_expected[148] => 0., 
xpopstr_expected[149] => 0, xpopstr_expected[150] => 0, 
xpopstr_expected[151] => -780.606, xpopstr_expected[152] => 0., 
xpopstr_expected[153] => 0, xpopstr_expected[154] => 0, 
xpopstr_expected[155] => 0., xpopstr_expected[156] => 1.21539, 
xpopstr_expected[157] => 0, xpopstr_expected[158] => 0, 
xpopstr_expected[159] => 0., xpopstr_expected[160] => -1.86789, 
xpopstr_expected[161] => 0, xpopstr_expected[162] => 0, 
xpopstr_expected[163] => 0., xpopstr_expected[164] => -24.2693, 
xpopstr_expected[165] => 0, xpopstr_expected[166] => 0, 
xpopstr_expected[167] => 0., xpopstr_expected[168] => -245.959, 
xpopstr_expected[169] => 0, xpopstr_expected[170] => 0, 
xpopstr_expected[171] => 0., xpopstr_expected[172] => -3199.36, 
xpopstr_expected[173] => 0, xpopstr_expected[174] => 0, 
xpopstr_expected[175] => 0.61636, xpopstr_expected[176] => 0, 
xpopstr_expected[177] => 0, xpopstr_expected[178] => 0, 
xpopstr_expected[179] => -3.94401, xpopstr_expected[180] => 0., 
xpopstr_expected[181] => 0, xpopstr_expected[182] => 0, 
xpopstr_expected[183] => -8.18325, xpopstr_expected[184] => 0., 
xpopstr_expected[185] => 0, xpopstr_expected[186] => 0, 
xpopstr_expected[187] => 19.0191, xpopstr_expected[188] => 0., 
xpopstr_expected[189] => 0, xpopstr_expected[190] => 0, 
xpopstr_expected[191] => 800.881, xpopstr_expected[192] => 0., 
xpopstr_expected[193] => 0, xpopstr_expected[194] => 0, 
xpopstr_expected[195] => 18280.6, xpopstr_expected[196] => 0., 
xpopstr_expected[197] => 0, xpopstr_expected[198] => 0, 
xpopstr_expected[199] => 0., xpopstr_expected[200] => 2.77362, 
xpopstr_expected[201] => 0, xpopstr_expected[202] => 0, 
xpopstr_expected[203] => 0., xpopstr_expected[204] => -9.48987, 
xpopstr_expected[205] => 0, xpopstr_expected[206] => 0, 
xpopstr_expected[207] => 0., xpopstr_expected[208] => -80.4034, 
xpopstr_expected[209] => 0, xpopstr_expected[210] => 0, 
xpopstr_expected[211] => 0., xpopstr_expected[212] => -623.083, 
xpopstr_expected[213] => 0, xpopstr_expected[214] => 0, 
xpopstr_expected[215] => 0., xpopstr_expected[216] => -5645.39, 
xpopstr_expected[217] => 0, xpopstr_expected[218] => 0, 
xpopstr_expected[219] => 1.34604, xpopstr_expected[220] => 0, 
xpopstr_expected[221] => 0, xpopstr_expected[222] => 0, 
xpopstr_expected[223] => -11.7121, xpopstr_expected[224] => 0., 
xpopstr_expected[225] => 0, xpopstr_expected[226] => 0, 
xpopstr_expected[227] => -3.89309, xpopstr_expected[228] => 0., 
xpopstr_expected[229] => 0, xpopstr_expected[230] => 0, 
xpopstr_expected[231] => 291.196, xpopstr_expected[232] => 0., 
xpopstr_expected[233] => 0, xpopstr_expected[234] => 0, 
xpopstr_expected[235] => 5240.32, xpopstr_expected[236] => 0., 
xpopstr_expected[237] => 0, xpopstr_expected[238] => 0, 
xpopstr_expected[239] => 95736.1, xpopstr_expected[240] => 0.
```
使用`reading-csv-nonlinear-SDP-x-power-11-standard-value.jl`读取`reading-csv-nonlinear-SDP-x-power-11-standard-value.jl`中的内容，保存至`nonlinear-SDP-x-power-11-standard-m-mat.jl`中。

将`nonlinear-SDP-x-power-11-standard-m-mat.jl`和最近一个block里面关于`xpopstr_expected`结合到`nonlinear-SDP-x-power-11-standard-point.jl`中。
但是命运和魔鬼不总是睡觉：结果当然是报错了：
```
ERROR: Feasibility checker for set type MathOptInterface.PositiveSemidefiniteConeTriangle has not been implemented yet.
```
无论如何，可以关闭SDP约束然后做检查。输出结果见`nonlinear-SDP-x-power-11-standard-point-check-oscillator-2-benchmark-with-ode.txt`。

能够观察到一些比较大的infeasible点，但是好像那个都是因为加减乘除过程中的舍入误差（是吗？）
最大的infeasible值对应如下两行：
```
-120 xpopstr_expected[76] - 600 xpopstr_expected[99] + 600 xpopstr_expected[124] + 200 xpopstr_expected[147] - 25 xpopstr_expected[172] - xpopstr_expected[195] + M[47,71] == 0.0    0.1738881602504989
-120 xpopstr_expected[76] - 600 xpopstr_expected[99] + 600 xpopstr_expected[124] + 200 xpopstr_expected[147] - 25 xpopstr_expected[172] - xpopstr_expected[195] + M[48,72] == 0.0    0.1738881602504989
```
对应的`M`分别是
```
M[47, 71] => -9782.922911839756,
M[48, 72] => -9782.922911839756
```
计算相对误差为`1.7774663238943776e-5`。这是否是浮点数计算误差？
另一个数据点是
```
-120 xpopstr_expected[99] + 240 xpopstr_expected[124] + 120 xpopstr_expected[147] - 20 xpopstr_expected[172] - xpopstr_expected[195] + M[60,60] == 0.0    0.1322533714774181
-120 xpopstr_expected[99] + 240 xpopstr_expected[124] + 120 xpopstr_expected[147] - 20 xpopstr_expected[172] - xpopstr_expected[195] + M[59,59] == 0.0    0.1322533714774181
```
对应的`M`是
```
M[59, 59] => -18013.035746628528,
M[60, 60] => -18013.035746628528
```
相对误差为`7.342092323453684e-6`。

我们来考虑一下`M[60, 60]`这个点。它对应的`M_index`是30, 30，于是运行以下诊断代码：
```julia
@show index_to_xpower(M_index_to_xpopstr_index[30])
@show index_to_ppower(M_index_to_xpopstr_index[30])
```
得到
```
index_to_xpower(M_index_to_xpopstr_index[30]) = 4
index_to_ppower(M_index_to_xpopstr_index[30]) = 5
```
唉其实也没必要检查什么东西。看着下面的式子：
```
-120 * (-18.9248) + 240 * -65.8857 + 120 * -117.932 - 20 * -3199.36 - 18280.6 -18013.035746628528
```
我就猜得出来怎么回事了。很明显xpopstr那些变量的有效数字保留少了……

## 2022.3.22

然后我们来诊断`M`矩阵：运行诊断代码
```julia
M_at_point = zeros(2 * (L_max + 1)^2, 2 * (L_max + 1)^2)
for i in 1 : 2 * (L_max + 1)^2
    for j in i : 2 * (L_max + 1)^2
        M_at_point[i, j] = point[M[i, j]]
        M_at_point[j, i] = point[M[i, j]]
    end
end

println(eigen(M_at_point).values)
```
然后会发现至少有两个问题
- 其中一个是，`calculate-all-correlation-functions-in-xp-oscillator-for-benchmark-3.nb`中
  $$
  \ii \to \begin{pmatrix}
    0 & -1 \\ 1 & 0
  \end{pmatrix}
  $$
  的做法是错误的
- 第二个问题更大，就是无论是我的Julia代码还是Mathematica代码都没有做$O_i^\dagger O_j$里面的转置
  恐怕这个才是真正大的问题……（这也解释了为什么只有$x$时正定性条件是成立的）

首先尝试修复Mathematica中的这个问题。
- 修复后，$M$矩阵的谱形如`M-matrix-benchmark-spectrum.pdf`
- 矩阵本身保存在`nonlinear-SDP-x-power-11-standard-m-value.csv`中（覆盖原文件）

修改
```julia
op_ij = xpopstr_normal_ord(op1_idx_xpower, op1_idx_ppower, op2_idx_xpower, op2_idx_ppower)
```
为
```julia
op_ij = xpopstr_normal_ord(0, op1_idx_ppower, op1_idx_xpower + op2_idx_xpower, op2_idx_ppower)
```
这么做了以后优化仍然失败，不过至少能够看到进展。

重复`reading-csv-nonlinear-SDP-x-power-11-standard-value.jl`中的步骤，覆盖`nonlinear-SDP-x-power-11-standard-m-mat.jl`和`nonlinear-SDP-x-power-11-standard-point.jl`。

重新做feasible测试。然后这回肯定是暴露出来一些问题了，因为有很多constraint是真的不满足。哦不对，是`nonlinear-SDP-x-power-11-standard-point.jl`中的测试部分写错了。

为了benchmark方便，将扩张后的$M$矩阵的特征值列举如下：
```
{152336. + 0. I, 152336. + 0. I, 35073.5 + 0. I, 35073.5 + 0. I, 
 7871.89 + 0. I, 7871.89 + 0. I, 3219.71 + 0. I, 3219.71 + 0. I, 
 1110.89 + 0. I, 1110.89 + 0. I, 418.283 + 0. I, 418.283 + 0. I, 
 143.947 + 0. I, 143.947 + 0. I, 142.996 + 0. I, 142.996 + 0. I, 
 34.4983 + 0. I, 34.4983 + 0. I, 25.0739 + 0. I, 25.0739 + 0. I, 
 10.7244 + 0. I, 10.7244 + 0. I, 1.07144 + 0. I, 1.07144 + 0. I, 
 0.186934 + 5.2322*10^-15 I, 0.186934 - 5.2322*10^-15 I, 
 0.0251751 + 9.67326*10^-14 I, 0.0251751 - 9.67326*10^-14 I, 
 0.0198789 + 0. I, 
 0.0198789 + 0. I, -0.0059278 + 0. I, -0.0059278 + 0. I, 
 0.00248624 + 0. I, 0.00248624 + 0. I, 0.0019754 + 0. I, 
 0.0019754 + 0. I, -0.00160649 + 0. I, -0.00160649 + 
  0. I, -0.00122542 + 0. I, -0.00122542 + 0. I, -0.001218 + 
  0. I, -0.001218 + 0. I, 0.000608623 + 0. I, 0.000608623 + 0. I, 
 0.000213875 + 1.46494*10^-14 I, 0.000213875 - 1.46494*10^-14 I, 
 0.000120688 + 0. I, 
 0.000120688 + 0. I, -0.0000893295 + 0. I, -0.0000893295 + 0. I, 
 0.0000463036 + 0. I, 
 0.0000463036 + 0. I, -0.0000309406 + 
  3.27116*10^-15 I, -0.0000309406 - 3.27116*10^-15 I, -0.0000208805 + 
  0. I, -0.0000208805 + 0. I, 
 1.053*10^-12 + 0. I, -7.07614*10^-14 + 0. I, -1.50663*10^-14 + 
  2.04839*10^-14 I, -1.50663*10^-14 - 2.04839*10^-14 I, 
 1.17326*10^-14 + 9.41112*10^-15 I, 
 1.17326*10^-14 - 9.41112*10^-15 I, -9.51811*10^-15 + 
  3.10301*10^-15 I, -9.51811*10^-15 - 3.10301*10^-15 I, 
 3.84313*10^-15 + 8.77461*10^-15 I, 3.84313*10^-15 - 8.77461*10^-15 I,
  2.66106*10^-15 + 6.21263*10^-15 I, 
 2.66106*10^-15 - 6.21263*10^-15 I, 
 5.67874*10^-15 + 0. I, -4.93221*10^-15 + 0. I, 
 3.09633*10^-15 + 0. I, -1.49068*10^-16 + 0. I}
```

在`nonlinear-SDP-x-power-11-standard-point.jl`中新增的两个测试
```julia
##

point[xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :real)] + 
    point[xpopstr_expected_real_imag_parts(xpopstr_index(0, 2), :real)] + 
    g * point[xpopstr_expected_real_imag_parts(xpopstr_index(4, 0), :real)]

##

norm(map(x -> point[x], M)' - map(x -> point[x], M))
```
均没有什么问题。因此现在的`jump-oscillator-2-benchmark-with-ode.jl`应该是正确的了？

于是当务之急是，看看它的最终优化结果是什么。

为了避免把我自己的电脑算爆炸，还是去超算上算比较好。
`print-test.jl`用来看在那个超算上输出是怎么记录的。

将`2022-3-22-run-1.sh`和`jump-oscillator-3.jl`提交，用前者控制计算。
核对：以上是全部需要提交的文件。`2022-3-22-run-1.sh`的是我需要计算的程序。
```sh
sc81233@ln71:~/jinyuanwu/numerical-bootstrap$ sbatch 2022-3-22-run-1.sh
Submitted batch job 880535
```
看输出应该是跑起来了。不过收敛得真的非常、非常慢。可能要将`max_iter`调到更大，，，

观察`slurm-880535.out`可以看到能量确实是一直在下降，但是收敛得非常，非常慢。

将`max_iter`调到一千万？

在`jump-oscillator-3.jl`中做对应的修改。

将`2022-3-22-run-2.sh`和`jump-oscillator-3.jl`提交，用前者控制计算。
核对：以上是全部需要提交的文件。`2022-3-22-run-1.sh`的是我需要计算的程序。
```sh
sc81233@ln71:~/jinyuanwu/numerical-bootstrap$ sbatch 2022-3-22-run-2.sh
Submitted batch job 880774
```

优化是优化完了，得到的数据是`slurm-880774.out`，但是最后几行是
```
3353475	 1.3968e+01	2.3154e+00	7.5731e-04	2.2861e-06
3353500	 1.3968e+01	2.0787e+00	1.9817e-05	2.2861e-06

------------------------------------------------------------------
>>> Results
Status: Solved
Iterations: 3353517 (incl. 17 safeguarding iter)
Optimal objective: 13.97
Runtime: 11580.482s (1.158048171e7ms)

-----------------------------------------------------------
Results:

Objective value:   13.96763393230745
x square expectation:     1.6702472405989577
```
这个很明显有问题啊，这个表的header是
```
Iter:	Objective:	Primal Res:	Dual Res:	Rho:
```

## 2022.3.23

我们有
```
julia> @show get_optimizer_attribute(model, "eps_abs")
get_optimizer_attribute(model, "eps_abs") = 1.0e-5
1.0e-5

julia> @show get_optimizer_attribute(model, "eps_rel")
get_optimizer_attribute(model, "eps_rel") = 1.0e-5
1.0e-5
```
结合COSMO的[收敛条件](https://oxfordcontrol.github.io/COSMO.jl/stable/method/#Termination-criteria)，引入以下配置：
```
set_optimizer_attributes(model, "max_iter" => 10000000, "eps_rel" => 1.0e-10)
```
保存在`jump-oscillator-3-cosmo.jl`里面。

将`jump-oscillator-3-cosmo.jl`和`2022-3-23-run-1.sh`复制到超算上，检查无误。然后等着吧，超算占用满了……

## 2022.3.24

将`jump-oscillator-3-cosmo.jl`改名为`jump-oscillator-3-cosmo-opsrel-1e-10.jl`，和`2022-3-24-run-1.pbs`复制到超算上，核对无误后运行

## 2022.3.25

将`2022-3-24-run-1.pbs`任务的结果保存于`2022-3-24-run-1-out`中。
顺带着把前几次的运行结果推到Github上面。out文件确实大了一点，是不是要用LFS？
```
remote: warning: See http://git.io/iEPt8g for more information.
remote: warning: File oscillator-simple-prototype/2022-3-24-run-1-out is 61.38 MB; this is larger than GitHub's recommended maximum file size of 50.00 MB
remote: warning: GH001: Large files detected. You may want to try Git Large File Storage - https://git-lfs.github.com.
To github.com:wujinq/numerical-boostrap.git
   b0b37b1..f1da28d  main -> main
```

不过`2022-3-24-run-1-out`确实体现出了一些比较糟糕的东西。优化能量到`4.7624e+00`附近以后，就长时间没有进展了。然后，能量开始缓慢上升。16933250迭代左右往后的计算——总计算时间的一小半——都完全浪费了。

目前可以继续做的事情包括：
- 更换求解器
- 固定几个参数，看看会有什么反应（无论如何，先做出来一个收敛的结果再说）
- 再次检查feasible条件

先从容易的做起吧。固定几个参数，看看有什么反应。

在`jump-oscillator-3-cosmo-opsrel-1e-10-fix-x2.jl`中固定$x^2$。使用`2022-3-25-run-1.pbs`提交此任务。
核对涉及到的文件
- `2022-3-25-run-1.pbs`
- `jump-oscillator-3-cosmo-opsrel-1e-10-fix-x2.jl`
  
无误。提交。

结果发现提示infeasible……难道还有隐藏的bug吗？？？文件见`2022-3-25-run-1.out`。

然后更换求解器。在`jump-oscillator-3-csdp.jl`中做这件事。使用使用`2022-3-25-run-2.pbs`提交此任务。核对无误后提交。
用CSDP发现直接报告lack of progress。

## 2022.3.26

那么，看起来我们需要再次检查是不是约束写错了。实际上就是再次检查feasibility。我们在`jump-oscillator-3-feasibility-checking.jl`里面做这件事。

运行诊断代码
```julia
value(xpopstr_expected_real_imag_parts(xpopstr_index(2, 0), :real) + 
    xpopstr_expected_real_imag_parts(xpopstr_index(0, 2), :real) + 
    g * xpopstr_expected_real_imag_parts(xpopstr_index(4, 0), :real), v -> point[v])
```
得到`1.3921990000000002`。因此能量没问题。

然后看`M`矩阵。运行诊断代码
```julia
eigen(value.(M, v -> point[v])).values
```
得到
```
72-element Array{Complex{Float64},1}:
  -0.005927795262277084 + 0.0im
 -0.0059277952622316015 + 0.0im
 -0.0016064864154372303 + 0.0im
 -0.0016064864154119894 + 0.0im
 -0.0012254218691151704 - 1.2406040499328298e-14im
 -0.0012254218691151704 + 1.2406040499328298e-14im
 -0.0012179953852214963 + 0.0im
 -0.0012179953851737479 + 0.0im
  -8.932946992171226e-5 + 0.0im
  -8.932946990707611e-5 + 0.0im
  -3.094063366678461e-5 - 1.4616356095333473e-14im
                        ⋮
     1110.8909248175623 + 0.0im
     1110.8909248175642 + 0.0im
     3219.7065089519874 + 0.0im
     3219.7065089519892 + 0.0im
      7871.886495860934 + 0.0im
      7871.886495860938 + 0.0im
      35073.54735146827 + 0.0im
      35073.54735146832 + 0.0im
      152336.1523407796 + 0.0im
     152336.15234077966 + 0.0im
```

在`calculate-all-correlation-functions-in-xp-oscillator-for-benchmark-4.nb`中获得更加精确的`xpopstr_expected`值。
根据结果，修改`nonlinear-SDP-x-power-11-standard-point.jl`。
这回，运行
```julia
max(collect(values(infeasible_report))...)
```
得到的就是`5.820766091346741e-11`，且infeasibility report的大小也缩水了很多。因此基本上可以确定是数值原因导致了之前的“大的infeasibility”。

运行
```julia
real_part_eigen_m = eigen(value.(M, v -> point[v])).values |> real
min(real_part_eigen_m[real_part_eigen_m .< 0]...)
```
得到-0.006左右的值。绝对值也不算小了，虽然所有特征值的乘积的绝对值非常非常小。与`calculate-all-correlation-functions-in-xp-oscillator-for-benchmark-4.nb`相比看不出明显区别。
因此，大概率我们又遇到了一个数值误差。

因此，基本可以确定，下面我们要做的事情是纯粹的优化技巧问题。

首先将`nonlinear-SDP-x-power-11-standard-point.jl`中的`point`分离出来，单独放在文件`nonlinear-SDP-x-power-11-standard-point-def.jl`里面。
看看能否固定变量，至少先跑出来一个feasible的解。

修改`jump-oscillator-3-cosmo-opsrel-1e-10-fix-x2.jl`。在`2022-3-26-run-1.sh`中调用它。

将如下文件：
- `jump-oscillator-3-cosmo-opsrel-1e-10-fix-x2.jl`
- `2022-3-26-run-1.sh`
- `nonlinear-SDP-x-power-11-standard-point-def.jl`

提交到服务器。核对无误。然后又是Batch job submission failed。匪夷所思。

更换服务器以后发现还是“infeasible”……