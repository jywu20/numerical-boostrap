开始做二维Hubbard model。

# 2022.5.10

发现即使$K=5$都会出现内存溢出。在建立COSMO模型时程序被系统杀掉了。

# 2022.5.12

目前的测试结果是，$K=8$时会出现“M矩阵里面的一些算符不在`hubbard_opstr_basis`中”的老错误。这个错误实际上说明`hubbard_opstr_basis`的筛选标准是有问题的，放过了很多本应该筛选进去的算符。
之前一维情况下，空间反演对称性竟然会弄出`hubbard_opstr_basis`以外的算符就很能说明问题了。

在`operator-algebra.jl`中加入了输出M矩阵内容的功能。我们发现在$K=8$时，`c(91) c(191)`没有出现在`hubbard_opstr_basis`中。
这真是稀奇，因为`c(101) c(191)`出现了！

为了便于分析，也许应该把所有的label也输出。

在`operator-algebra.jl`中，`c(91) c(191)`从一开始就没有出现在`qualified_opstr_annihilate_half`中。

我真的是弱智，真的是。`check_max_iter_opstr_num`不关。

# 2022.5.13

现在要做的是至少跑出一些结果来，好先把文章填满。

然后还要施加旋转对称性。