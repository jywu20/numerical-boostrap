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