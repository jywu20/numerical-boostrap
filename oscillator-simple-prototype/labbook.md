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

算了我还是老老实实解析算吧……
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