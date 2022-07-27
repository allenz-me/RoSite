---
title: "第四章 分布鲁棒优化（模糊集、机会约束问题、分布鲁棒线性优化）"
date: 2022-01-06
# menu: "main"
draft: false
# you can close something for this content if you open it in config.toml.
mathjax: false  # 手动加入了 mathjax js 脚本
toc: true
---



现实世界的优化问题指在满足相关约束条件的前提下，确定一组决策变量的值，使预设的目标函数值最优。相关研究成果（理论、模型、算法、应用）在管理科学、金融工程、军事指挥等领域发挥着巨大指导作用，创造了巨大的社会经济价值。

由于本书定位为入门级、科普级，本章将以优化问题中一类简单但重要的线性优化（亦称线性规划）问题为例，阐述**分布鲁棒优化**的基本思想、基本模型、基本结论。感兴趣的读者需细读相关文献以获得更深入广泛的理解。 

线性规划问题 (4.1) 中，决策变量为 $\boldsymbol{x} \in \mathbb{R}^I$，环境参数包括费用向量 $\boldsymbol{a}_0 \in \mathbb{R}^I$、约束条件左端项系数向量 $\boldsymbol{a}_m \in \mathbb{R}^I$、右端项系数 $b_m \in \mathbb{R}$，这些参数为确定值。线性规划可用于解决许多现实问题，例如投资组合优化、生产计划、最短路等问题。
$$
\begin{equation}
\label{ro.mod.lo}
\begin{array}{rll}
\displaystyle \min_{\boldsymbol{x}} & \boldsymbol{a}_0^\top \boldsymbol{x}, \\
\mbox{s.t.} & \boldsymbol{a}_m^\top \boldsymbol{x} \le b_m, & m \in [M]. \tag{4.1}
\end{array}
\end{equation}
$$
现实世界中描述未来发生事件的环境参数在优化/规划/计划阶段往往不确定，例如未来某商品需求量、两地间旅行时长、某股票回报率等。为了让优化结果对现实更具指导意义，在优化模型中考虑环境参数的不确定性至关重要。

在不确定环境下，(4.1) 中对于某一 $m \in [M]$ 的约束式变成了
$$
\begin{equation}
\label{dro.con.1}
\boldsymbol{a}(\tilde{\boldsymbol{\varepsilon}})^\top\boldsymbol{x} \le  b(\tilde{\boldsymbol{\varepsilon}}). \tag{4.2}
\end{equation}
$$
其中，为了阐述方便，忽略下标 $m$；随机变量 $\tilde{\boldsymbol{\varepsilon}}$ 表示影响环境参数的随机因素（例如，旅行时长受天气、交通灯时长等随机因素影响），假设 $\boldsymbol{a}(\tilde{\boldsymbol{\varepsilon}})$ 和 $b(\tilde{\boldsymbol{\varepsilon}})$ 皆为 $\tilde{\boldsymbol{\varepsilon}}$ 的仿射函数，即
$
\boldsymbol{a}(\tilde{\boldsymbol{\varepsilon}}) \triangleq \boldsymbol{a}^0 + \sum_{j \in [J]}\boldsymbol{a}^j \tilde{\varepsilon}_j, b(\tilde{\boldsymbol{\varepsilon}}) \triangleq b^0 + \sum_{j \in [J]}b^j \tilde{\varepsilon}_j,
$
则有
$$
\begin{equation*}
\boldsymbol{a}(\tilde{\boldsymbol{\varepsilon}})^\top \boldsymbol{x} - b(\tilde{\boldsymbol{\varepsilon}})
= \underbrace{(\boldsymbol{a}^{0})^\top\boldsymbol{x} - b^0}_{=y^0(\boldsymbol{x})} + \sum_{j \in [J]}\big(\underbrace{(\boldsymbol{a}^{j})^\top\boldsymbol{x} - b^j}_{=y^j(\boldsymbol{x})}\big)\tilde{\varepsilon}_j = y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})'\tilde{\boldsymbol{\varepsilon}}.
\end{equation*}
$$
因此，约束式 (4.2) 等价于
$$
\begin{equation}
\label{dro.con.2}
y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top\tilde{\boldsymbol{\varepsilon}} \le 0. \tag{4.3}
\end{equation}
$$
不确定环境下的线性规划从技术上主要关注如何处理约束式 (4.3)。因其左端项为随机变量而右端项为实数，故通常意义上无法直接比较大小，“$\le$” 符号用在此处不够严谨。对此，本章主要介绍两种典型处理方式，第4.2节基于**分布鲁棒机会约束规划**思想，介绍如何处理
$$
\begin{equation}
\label{dro.con.cc}
\mathbb{P} [y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top\tilde{\boldsymbol{\varepsilon}} \le 0] \ge 1 - \epsilon, \quad \forall \mathbb{P} \in \mathscr{F}, \tag{4.4}
\end{equation}
$$
也就是约束 (4.3) 成立的概率不小于 $1 - \epsilon$，其中阈值 $\epsilon \in [0, 1]$ 典型取值为 1\% 或 5\%。第4.3节基于**分布鲁棒线性优化**范式，介绍如何处理
$$
\begin{equation}
\label{dro.con.lo}
\mathbb{E}_{\mathbb{P}} [y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top\tilde{\boldsymbol{\varepsilon}}] \le 0, \quad \forall \mathbb{P} \in \mathscr{F}, \tag{4.5}
\end{equation}
$$

也就是约束 (4.3) 需在其左端项通过均值来度量的情况下满足。

事实上，目标函数参数不确定性亦可纳入约束讨论，因为通过引入辅助决策变量 $t \in \mathbb{R}$，(4.1)等价于
$$
\begin{array}{rll}
\displaystyle \min_{t, \boldsymbol{x}} & t, \\
\mbox{s.t.} & \boldsymbol{a}_0^\top \boldsymbol{x} \le t, \\
&\boldsymbol{a}_m^\top \boldsymbol{x} \le b_m, & m \in [M]. 
\end{array}
$$
进而将目标函数中不确定参数 $\boldsymbol{a}_0$ 置于约束中。

约束式 (4.4) 和 (4.5) 中，随机变量 $\tilde{\boldsymbol{\varepsilon}}$ 服从联合概率分布 $\mathbb{P}$，而 $\mathbb{P}$ 本身也不确定，属于模糊集 $\mathscr{F}$；第4.1节将介绍模糊集相关内容。分布鲁棒优化采取保守策略，令这两个约束条件对模糊集中所有概率分布皆满足，它也是因此得名---“鲁棒”的内涵是考虑最坏情况，而“分布”表明最坏情况的主体是环境参数的分布函数。望本章内容能抛砖引玉，启发读者研究和处理更复杂的形式，并用于解决实际问题。

## 4.1  模糊集（Ambiguity set）

问题环境中随机参数的分布函数往往难以从现实世界直接获取。鉴于此，分布鲁棒优化方法假设其分布函数并不明确，而是处于一个**模糊集**(ambiguity set) 中。模糊集通过随机变量的不完全分布信息构建而成。特别地，它还需保证相应分布鲁棒优化模型在计算上可处理 (tractable)，也就是现实规模问题可在允许时间范围内求解。

从数学上说，分布鲁棒优化囊括随机规划(stochastic programming) 和传统鲁棒优化(robust optimization)为特殊形式，因此分布鲁棒优化更具一般性。当环境变量 $\tilde{\boldsymbol{\varepsilon}}$ 的分布函数 $\mathbb{P}_0$ 可获知时，可令模糊集为单元素集 $\mathscr{F}_S \triangleq \{\mathbb{P}_0\}$，则分布鲁棒优化退化为随机规划；当仅知环境变量的不确定集 $\Xi$ 时，可令模糊集为 $\mathscr{F}_R \triangleq \mathscr{P}_0(\Xi)$，即支撑集为 $\Xi$ 的所有概率分布函数之集合，则分布鲁棒优化退化为经典鲁棒优化。

按照描述分布函数的信息种类划分，目前相关研究主要提出了两类模糊集。

###　4.1.1  基于广义矩信息（generalized moment information）的模糊集

在统计学中，矩 (moment) 表征随机变量的分布。对于随机变量 $\tilde{\varepsilon}$，其 $n$ 阶矩被定义为 $\mathbb{E}_\mathbb{P}[\tilde{\varepsilon}^n]$，$n \ge 1$。因此，随机变量一阶矩为均值，表征其位置 (location)，二阶矩与方差有关，表征其散度 (dispersion)，三阶矩表征其偏斜度，等等。

更广义地，还可利用其它形式表征随机变量的位置、散度、偏斜度等特性。例如，绝对离差均值  $\mathbb{E}_{\mathbb{P}}[|\tilde{\varepsilon} - \mu|]$ 可表征 $\tilde{\varepsilon}$ 的散度，其中 $\mu$ 为其均值。再如，半绝对离差均值 $\mathbb{E}_{\mathbb{P}}[(\tilde{\varepsilon} - \mu)^+]$ 和 $\mathbb{E}_{\mathbb{P}}[(\mu - \tilde{\varepsilon})^+]$ 可从某种程度刻画 $\tilde{\varepsilon}$ 的偏斜度，其中 $(x)^+ \triangleq \max\{x, 0\}$。

早期研究往往假设随机参数的概率分布无法准确获取，但其部分广义矩信息（和支撑集）可获取或估计，于是根据这些信息构建模糊集。例如，通过 $\tilde{\boldsymbol{\varepsilon}}$ 的均值 $\boldsymbol{\mu}$ 和协方差矩阵 $\boldsymbol{\Sigma}$ 构成的模糊集 (Ghaoui et al., 2003; Popescu, 2007; Chen and Sim, 2009)  为
$$
\begin{equation}
\label{eq.dro.as.mv}
\mathscr{F}_{MV} = \left\{
\mathbb{P} \in \mathscr{P}_0(\mathbb{R}^J)
\left|\begin{array}{l}
\tilde{\boldsymbol{\varepsilon}} \sim \mathbb{P} \\
\mathbb{E}_\mathbb{P}[\tilde{\boldsymbol{\varepsilon}}] = \boldsymbol{\mu} \\
\mathbb{E}_\mathbb{P}[(\tilde{\boldsymbol{\varepsilon}} -\boldsymbol{\mu}) (\tilde{\boldsymbol{\varepsilon}} - \boldsymbol{\mu})'] = \boldsymbol{\Sigma} \\ \tag{4.6}
\end{array}
\right.
\right\}.
\end{equation}
$$
如果进一步考虑支撑集 $\Xi$，则模糊集为 $\mathscr{F}_{MVS} = \mathscr{F}_{MV} \cap \mathscr{P}_0(\Xi)$. 但研究表明，基于 $\mathscr{F}_{MVS}$ 的分布鲁棒优化模型一般不可处理 (Bertsimas and Popescu, 2005; Natarajan et al., 2011)。而如果给定的 $\boldsymbol{\Sigma}$ 不是准确协方差而是协方差的上界时，则模糊集为
$$
\begin{equation}
\label{eq.dro.as.moment}
\mathscr{F}_{M} = \left\{
\mathbb{P} \in \mathscr{P}_0(\Xi)
\left|
\begin{array}{l}
\tilde{\boldsymbol{\varepsilon}} \sim \mathbb{P} \\
\mathbb{E}_\mathbb{P}[\tilde{\boldsymbol{\varepsilon}}] = \boldsymbol{\mu} \\
\mathbb{E}_\mathbb{P}[(\tilde{\boldsymbol{\varepsilon}} -\boldsymbol{\mu}) (\tilde{\boldsymbol{\varepsilon}} - \boldsymbol{\mu})'] \preceq \boldsymbol{\Sigma}
\end{array} \tag{4.7}
\right.
\right\},
\end{equation}
$$
其中 “$\preceq$”为半正定锥空间意义上的小于等于，也就是 $\boldsymbol{X} \preceq \boldsymbol{Y}$ 意味着 $\boldsymbol{Y} - \boldsymbol{X}$ 为半正定矩阵。有趣的是，基于 $\mathscr{F}_{M}$ 的分布鲁棒线性优化模型却可处理(Wiesemann et al., 2014; Hanasusanto et al., 2015)。

Delage and Ye (2010) 研究了 $\mathscr{F}_{M}$ 的一个变种
$$
\begin{equation}
\label{eq.dro.as.moment.dy}
\mathscr{F}_{DY} = \left\{
\mathbb{P} \in \mathscr{P}_0(\Xi)
\left|
\begin{array}{l}
\tilde{\boldsymbol{\varepsilon}} \sim \mathbb{P} \\
(\mathbb{E}_{\mathbb{P}}[\tilde{\boldsymbol{\varepsilon}}] - \boldsymbol{\mu})^\top \boldsymbol{\Sigma}^{-1} (\mathbb{E}_{\mathbb{P}}[\tilde{\boldsymbol{\varepsilon}}] - \boldsymbol{\mu}) \le \gamma_1 \\
\mathbb{E}_\mathbb{P}[(\tilde{\boldsymbol{\varepsilon}} -\boldsymbol{\mu}) (\tilde{\boldsymbol{\varepsilon}} - \boldsymbol{\mu})'] \preceq \gamma_2 \boldsymbol{\Sigma}
\end{array} \tag{4.8}
\right.
\right\},
\end{equation}
$$
其中第一个约束指 $\tilde{\boldsymbol{\varepsilon}}$ 的均值处于一个以 $\boldsymbol{\mu}$ 为球心的椭球中，$\gamma_1 \ge 0$ 和 $\gamma_2 \ge 1$ 为两个参数。从数据驱动的视角看，假设 $\tilde{\boldsymbol{\varepsilon}}$ 客观上服从概率分布 $\mathbb{P}_0$，但无法观测该分布，而仅能观测其 $N$ 组样本/历史数据/观测值 $(\hat{\boldsymbol{\varepsilon}}_\omega)_{\omega \in [N]}$，令 
$$
\boldsymbol{\mu} \triangleq \frac{1}{N}\sum_{\omega \in [N]} \hat{\boldsymbol{\varepsilon}}_\omega, \qquad \boldsymbol{\Sigma} \triangleq \frac{1}{N} \sum_{\omega \in [N]} (\hat{\boldsymbol{\varepsilon}}_\omega - \boldsymbol{\mu}) (\hat{\boldsymbol{\varepsilon}}_\omega - \boldsymbol{\mu})^\top,
$$
且 $\gamma_1$ 和 $\gamma_2$ 通过与样本量 $N$ 和参数 $\delta > 0$ 有关的某函数给定时（随着 $N \rightarrow \infty$，有 $\gamma_1 \rightarrow 0$ 和 $\gamma_2 \rightarrow 1$），则 Delage and Ye (2010) 证明了在统计学上 $\mathbb{P}_0 \in \mathscr{F}_{DY}$ 的置信度大于等于 $1 - \delta$。

并非任意基于广义矩信息（和支撑集）的模糊集都能保证相应分布鲁棒优化模型可处理。Wisemann，Kuhn 和 Sim 提出了一种具有一般性的模糊集表达形式，能囊括 $\mathscr{F}_M$ 和 $\mathscr{F}_{DY}$ 为其特殊形式，能建模许多其它的广义矩信息，如绝对离差、半方差、高阶矩等，且（在一些技术性假设条件下）对于分布鲁棒线性优化在计算上可处理 (Wiesemann et al., 2014; Hanasusanto et al., 2015)。其形式为
$$
\begin{equation}
\label{eq.dro.as.wks}
\mathscr{F}_{WKS} = \left\{
\mathbb{P} \in \mathscr{P}_0(\mathbb{R}^J \times \mathbb{R}^L)
\left|
\begin{array}{ll}
(\tilde{\boldsymbol{\varepsilon}}, \tilde{\boldsymbol{u}}) \sim \mathbb{P} \\
\mathbb{E}_\mathbb{P}[\boldsymbol{A}\tilde{\boldsymbol{\varepsilon}} + \boldsymbol{B} \tilde{\boldsymbol{u}}] = \boldsymbol{b} \\
\mathbb{P}[(\tilde{\boldsymbol{\varepsilon}}, \tilde{\boldsymbol{u}}) \in \Xi_k] \in [\underline{p}_k, \overline{p}_k], \forall k \in [K]
\end{array}
\right.
\right\},
\end{equation} \tag{4.9}
$$
其中 $\mathbb{P}$ 为 $\tilde{\boldsymbol{\varepsilon}}$ 和辅助随机变量 $\tilde{\boldsymbol{u}}$ 的联合概率分布，$\boldsymbol{A} \in \mathbb{R}^{J \times Q}$，$\boldsymbol{B} \in \mathbb{R}^{L \times Q}$，$\boldsymbol{b} \in \mathbb{R}^Q$；置信集合 $\Xi_k$ 给定为
$$
\begin{equation}
\label{eq.support.set}
\Xi_k = \{
(\boldsymbol{\varepsilon}, \boldsymbol{u}) \in \mathbb{R}^J \times \mathbb{R}^L |
\boldsymbol{C}_k \boldsymbol{\varepsilon} + \boldsymbol{D}_k \boldsymbol{u} \preceq_{\mathscr{K}_k} \boldsymbol{c}_k
\},
\end{equation} \tag{4.10}
$$

其中 $\boldsymbol{C}_k \in \mathbb{R}^{J \times R}$，$\boldsymbol{D}_k \in \mathbb{R}^{L \times R}$，$\boldsymbol{c} \in \mathbb{R}^R$，而 $\mathscr{K}_k$ 代表某一真锥 (proper cone)，如非负象限、二阶锥、半正定锥等。在此，“$\preceq_{\mathscr{K}_k}$”是在该锥空间意义上的小于等于；关于锥和相应的锥规划 (conic programming) 问题相关介绍，可参见  Ben-Tal and Nemirovski (2001)。对于表示概率界的 $\underline{\boldsymbol{p}},\overline{\boldsymbol{p}} \in [0, 1]^K$，有 $\underline{\boldsymbol{p}} \le \overline{\boldsymbol{p}}$。

模糊集 $\mathscr{F}_{WKS}$ 巧妙之处在于引入了辅助随机变量 $\tilde{\boldsymbol{u}}$，这为计算上的可处理性提供了一种有效途径，如下例所示。

**示例4.1**：通过 $\tilde{\boldsymbol{\varepsilon}}$ 的均值 $\boldsymbol{\mu}$ 、绝对离差均值上界 $\boldsymbol{\sigma}$、半绝对离差均值上界 $\boldsymbol{h}$ 及满足 (4.10)形式的支撑集 $\Xi$ 构成的模糊集为
$$
\begin{equation}
\label{eq.dro.as.mad}
\mathscr{F}_{MAD} = \left\{
\mathbb{P} \in \mathscr{P}_0(\Xi)
\left|
\begin{array}{l}
\tilde{\boldsymbol{\varepsilon}} \sim \mathbb{P} \\
\mathbb{E}_\mathbb{P}[\tilde{\boldsymbol{\varepsilon}}] = \boldsymbol{\mu} \\
\mathbb{E}_\mathbb{P}[|\tilde{\boldsymbol{\varepsilon}} -\boldsymbol{\mu}|] \le \boldsymbol{\sigma} \\
\mathbb{E}_\mathbb{P}[(\tilde{\boldsymbol{\varepsilon}} -\boldsymbol{\mu})^+] \le \boldsymbol{h}
\end{array}
\right.
\right\},
\end{equation} \tag{4.11}
$$
其中 $|\tilde{\boldsymbol{\varepsilon}}|$ 代表对向量 $\tilde{\boldsymbol{\varepsilon}}$ 的每个元素分别取绝对值构成的向量，取正符 $(\tilde{\boldsymbol{\varepsilon}})^+$ 亦然。目前虽无法直接处理基于 $\mathscr{F}_{MAD}$ 的分布鲁棒优化模型，但如引入辅助变量，$\mathscr{F}_{MAD}$ 则变成
$$
\begin{equation}
\label{eq.dro.as.ma}
\mathscr{F}_{MADL} = \left\{
\mathbb{P} \in \mathscr{P}_0(\bar{\Xi})
\left|
\begin{array}{l}
(\tilde{\boldsymbol{\varepsilon}}, \tilde{\boldsymbol{u}}, \tilde{\boldsymbol{v}}) \sim \mathbb{P} \\
\mathbb{E}_\mathbb{P}[\tilde{\boldsymbol{\varepsilon}}] = \boldsymbol{\mu} \\
\mathbb{E}_\mathbb{P}[\tilde{\boldsymbol{u}}] = \boldsymbol{\sigma} \\
\mathbb{E}_\mathbb{P}[\tilde{\boldsymbol{v}}] = \boldsymbol{h}
\end{array}
\right.
\right\},
\end{equation} \tag{4.12}
$$

其中，扩展的支撑集为
$$
\bar{\Xi} = \{(\boldsymbol{\varepsilon}, \boldsymbol{u}, \boldsymbol{v}) | \boldsymbol{\varepsilon} \in \Xi,  \boldsymbol{u} \ge \boldsymbol{\varepsilon} -\boldsymbol{\mu}, \boldsymbol{u} \ge \boldsymbol{\mu} - \boldsymbol{\varepsilon}, \boldsymbol{v} \ge \boldsymbol{\varepsilon} - \boldsymbol{\mu}, \boldsymbol{v} \ge \boldsymbol{0}\}.
$$
在此，$\mathscr{F}_{MADL}$ 即为 $\mathscr{F}_{WKS}$ 的一个特例，因此变得可处理；对于具体处理方法，可参见 Wiesemann et al. (2014); Hanasusanto et al. (2015)。类似地，前文提到的 $\mathscr{F}_{DY}$ 和 $\mathscr{F}_{M}$ 以及文献中更多有趣的形式 （例如，Bertsimas et al., 2019），也可借助辅助变量化为 $\mathscr{F}_{WKS}$ 的形式。

### 4.1.2  基于统计距离 (statistical distance) 的模糊集

在大数据时代背景下，问题环境中随机参数的历史数据越来越容易获取。为了获得随机参数概率分布，一种自然的想法是通过历史数据对应的经验分布来近似描述真实概率分布。经验分布是建立在 $N$ 条历史数据点上的离散均匀分布，它视每一条历史数据 $\hat{\boldsymbol{\varepsilon}}_\omega$ 为随机变量的一个支撑点，出现概率为$1/N$，即
$$
\hat{\mathbb{P}}[\tilde{\boldsymbol{\varepsilon}}^\dagger = \hat{\boldsymbol{\varepsilon}}_\omega] = \frac{1}{N}, ~ \forall \omega \in [N].
$$
其中 $\tilde{\boldsymbol{\varepsilon}}^\dagger$ 表示 $\tilde{\boldsymbol{\varepsilon}}$ 对应的经验随机变量。

但经验分布并不等同于真正概率分布，分布鲁棒优化的思想是假设真正概率分布与经验分布在概率空间中的统计距离不超过某一阈值，并以此构建模糊集。
基于此思想，如何定义两个概率分布间的统计距离成为了关键，它不仅需要具有良好的统计学意义，而且需要保证相应的分布鲁棒优化模型可处理。在此例举两种典型形式。

第一是基于$\phi$-散度 ($\phi$-divergence) 的模糊集，定义为
$$
\begin{equation}
\label{eq.dro.as.pd}
\mathscr{F}_{\phi} = \left\{
\mathbb{P} \in \mathscr{P}_0(\mathbb{R}^J)
\left|
\begin{array}{l}
\tilde{\boldsymbol{\varepsilon}} \sim \mathbb{P}, \tilde{\boldsymbol{\varepsilon}}^\dagger \sim \hat{\mathbb{P}}, \\
D_\phi(\mathbb{P}||\hat{\mathbb{P}}) \le \theta
\end{array}
\right.
\right\},
\end{equation} \tag{4.13}
$$
其中 $D_\phi(\mathbb{P}||\hat{\mathbb{P}})$ 表示所认为的真正分布 $\mathbb{P}$ 对于经验分布 $\hat{\mathbb{P}}$ 的 **$\phi$-散度**（“距离”），定义如 (4.14)；而 $\theta > 0$ 为给定的“距离”上界。

$$
\begin{equation}
D_\phi(\mathbb{P}||\hat{\mathbb{P}}) = \sum_{\omega \in [N]} \hat{\mathbb{P}}(\omega) \phi \left(\frac{\mathbb{P}(\omega)}{\hat{\mathbb{P}}(\omega)}\right) \tag{4.14}
\end{equation}
$$

在 (4.14)中，$\phi:\mathbb{R}_+ \mapsto \mathbb{R}$ 为满足以下条件的凸函数：$\phi(1) = 0$，对于 $x > 0$ 有 $0\phi(x/0)\triangleq x \lim_{t \rightarrow +\infty} \phi(t) / t$，且 $0\phi(0 / 0) \triangleq 0$；$\mathbb{P}(\omega)$ 表示概率分布 $\mathbb{P}$ 中第 $\omega \in [N]$ 个观测值发生的概率。可见，该模糊集要求真正的概率分布函数支撑集与经验分布支撑集相同，也就无法考量到历史数据以外的点/场景，而仅仅是将历史数据发生的概率从经验分布的各 $1/N$ 变成了更具“鲁棒性”的值。

示例4.2：在此例举一种研究相对较多的 $\phi$ 函数形式：当 $\phi(x) = x \log x - x + 1$ 时，$\phi$散度具体化为 (4.15)，名为**Kullback-Leibler 散度**(KL divergence)，又名**相对熵 (relative entropy)**。
$$
\begin{equation}
D_{KL}(\mathbb{P}||\hat{\mathbb{P}}) = \sum_{\omega \in [N]} \mathbb{P}(\omega) \log \left(\frac{\mathbb{P}(\omega)}{\hat{\mathbb{P}}(\omega)}\right) \tag{4.15}
\end{equation}
$$
为后文阐述方便，将其模糊集记为 $\mathscr{F}_{KL}$，也就是 (4.13)中的将 $D_\phi$ 具体化为 $D_{KL}$。对于更多的 $\phi$ 函数形式，可参见 Ben-Tal et al. (2013); Jiang and Guan (2016)。

第二是基于 Wasserstein 距离的模糊集，定义为
$$
\begin{equation}
\mathscr{F}_W = \left\{ \mathbb{P} \in \mathscr{P}_0(\Xi)
~\left|~
\begin{array}{l}
\displaystyle \tilde{\boldsymbol{\varepsilon}} \sim \mathbb{P}, \tilde{\boldsymbol{\varepsilon}}^\dagger \sim \hat{\mathbb{P}}, \\
\displaystyle d_W(\mathbb{P}, \hat{\mathbb{P}}) \le \theta, \\
\end{array}\right.
\right\}. \tag{4.16}
\end{equation}
$$
此模糊集囊括了概率空间中以 Wasserstein 距离为度量标准，以经验分布 $\hat{\mathbb{P}}$ 为球心，以 $\theta \in \mathbb{R}_+$ 为半径的球中所有的概率分布。Esfahani and Kuhn (2018)的研究表明，记真实但未知的概率分布为 $\mathbb{P}_0$，则当 $\theta$ 通过与样本量 $N$ 和参数 $\mathbb{E}ta \in (0,1)$ 有关的某函数取值时（随着 $N \rightarrow \infty$，有 $\theta \rightarrow 0$），从统计学上可证明 $\mathbb{P} \in \mathscr{F}_W$ 的置信度大于等于 $1 - \mathbb{E}ta$。**Wasserstein 距离** $d_W:\mathcal{P}_0(\Xi) \times \mathcal{P}_0(\Xi) \mapsto [0, +\infty)$ 表示所考虑的分布与经验分布在概率空间中的一种距离，定义为
$$
\begin{equation}
\begin{array}{rll}
d_W(\mathbb{P}, \hat{\mathbb{P}})  = \inf & \displaystyle \mathbb{E}_{\bar{\mathbb{P}}} \big[\lVert \tilde{\boldsymbol{\varepsilon}} - \tilde{\boldsymbol{\varepsilon}}^\dagger \rVert \big] \\
\mbox{s.t.}	& \displaystyle \big(\tilde{\boldsymbol{\varepsilon}}, \tilde{\boldsymbol{\varepsilon}}^\dagger\big) \sim \bar{\mathbb{P}}, \\
& \displaystyle \tilde{\boldsymbol{\varepsilon}} \sim \mathbb{P}, \\
& \displaystyle \tilde{\boldsymbol{\varepsilon}}^\dagger \sim \hat{\mathbb{P}}, \\
& \displaystyle \bar{\mathbb{P}}\big[ (\tilde{\boldsymbol{\varepsilon}}, \tilde{\boldsymbol{\varepsilon}}^\dagger) \in \Xi \times \Xi \big] = 1,
\end{array}
\end{equation} \tag{4.17}
$$
其中的 $\bar{\mathbb{P}}$ 表示 $\tilde{\boldsymbol{\varepsilon}}$ 和 $\tilde{\boldsymbol{\varepsilon}}^\dagger$ 的\emph{联合概率分布}，$\lVert \cdot \rVert$ 表示范数。根据定义，可直观地将 Wasserstein 距离视为从真实分布 $\mathbb{P}$ 向经验分布 $\hat{\mathbb{P}}$ 移动概率质量 (probability mass) 的最小费用。上述定义准确说是**1型 Wasserstein 距离**，对于更一般的 Wasserstein 距离定义及更详细深入的介绍可参见Esfahani and Kuhn (2018); Gao and Kleywegt (2016); Zhao and Guan (2018)。

## 4.2 机会约束问题（Chance constraint）

机会约束规划是指当优化问题环境参数为随机变量时，在以一定概率满足约束条件的情况下进行优化。自从Charnes and Cooper (1959) 提出来以来，该框架在管理科学等各领域得到了广泛研究与应用。作为抛砖引玉，本节讨论如何处理分布鲁棒、独立机会约束(4.4)，读者可进而自行研究更一般也更难处理的联合机会约束 (4.18)：
$$
\begin{equation}
\mathbb{P} [y_m^0(\boldsymbol{x}) + \boldsymbol{y}_m(\boldsymbol{x})^\top\tilde{\boldsymbol{\varepsilon}} \le 0, ~\forall m \in [M]] \ge 1 - \epsilon, \quad \forall \mathbb{P} \in \mathscr{F},
\end{equation} \tag{4.18}
$$

其中的 $M$ 个约束同时成立的概率不小于 $1 - \epsilon$。

从计算角度看，约束式 (4.4) 左端项可等价表示为
$$
\begin{equation}
\mathbb{E}_\mathbb{P}\big[\mathbb{1}\{y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top\tilde{\boldsymbol{\varepsilon}} \le 0\}\big],
\end{equation} \tag{4.19}
$$
其中，$\mathbb{1}\{\cdot\}$ 为指示函数，当其事件发生取值为1，否则为0。指示函数为非凸函数，导致 (4.19) 对 $\boldsymbol{x}$ 或 $\tilde{\boldsymbol{\varepsilon}}$ 而言皆为非凸函数，除个别特殊情况（如 $\tilde{\boldsymbol{\varepsilon}}$ 服从联合正态分布）外难以处理。事实上，即便给定决策变量 $\boldsymbol{x}$ 和概率分布 $\mathbb{P}$，计算 (4.4) 左端项的概率值一般而言已是 NP 难问题，更何况还要基于此对 $\boldsymbol{x}$ 进行优化 (Nemirovski and Shapiro, 2006)。而有趣的是，在给定某些模糊集 $\mathscr{F}$ 的情况下，(4.4) 却可处理。

接下来讲述如何在分布鲁棒优化框架下对机会约束式 (4.4) 进行处理。易知，约束式 (4.4) 等价于
$$
\begin{equation}
\displaystyle \inf_{\mathbb{P} \in \mathscr{F}}\mathbb{P}[y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top\tilde{\boldsymbol{\varepsilon}} \le 0] \ge 1 - \epsilon. \tag{4.20}
\end{equation}
$$
当 (4.20) 中模糊集取 $\mathscr{F} \triangleq \mathscr{F}_{MV}$ 且其中的 $\boldsymbol{\mu} \triangleq \boldsymbol{0}$ 时， Ghaoui et al. (2003) 的研究表明，(4.20) 等价于
$$
\begin{equation}
\label{eq.dro.cc.mv.eq}
\displaystyle y^0(\boldsymbol{x}) + \sqrt{\frac{1 - \epsilon}{\epsilon}} \sqrt{\boldsymbol{y}(\boldsymbol{x})^\top\boldsymbol{\Sigma} \boldsymbol{y}(\boldsymbol{x})} \le 0. \tag{4.21}
\end{equation}
$$
这里，令 $\boldsymbol{\mu} \triangleq \boldsymbol{0}$ 并不失一般性，因为如果 $\boldsymbol{\mu} \not= \boldsymbol{0}$ 则可通过变量替换的方法，令 $\tilde{\boldsymbol{\xi}} \triangleq \tilde{\boldsymbol{\varepsilon}} - \boldsymbol{\mu}$ 并将 $\tilde{\boldsymbol{\xi}}$ 视为 (4.20) 中的 $\tilde{\boldsymbol{\varepsilon}}$。有趣的是，(4.21) 恰好等价于传统鲁棒优化约束式
$$
\begin{equation}
\label{eq.ro}
\displaystyle y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top \boldsymbol{\varepsilon} \le 0, \quad \forall \boldsymbol{\varepsilon} \in \Xi(\epsilon), \tag{4.22}
\end{equation}
$$

其中，不确定集为椭球形，给定为
$$
\Xi(\epsilon) \triangleq \left\{\boldsymbol{\varepsilon} \in \mathbb{R}^J \left|~\lVert \boldsymbol{\Sigma}^{1/2}\boldsymbol{\varepsilon}\rVert_2 \le \sqrt{\frac{1-\epsilon}{\epsilon}} \right. \right\}.
$$
关于 (4.21)与(4.22)的关系，可参见 Natarajan et al. (2009)。

当 (4.20) 中模糊集取 $\mathscr{F} \triangleq \mathscr{F}_{MVS}$ 时，分布鲁棒机会约束规划一般不可处理。研究者们提出了基于**条件风险值** (Conditional Value-at-Risk) 的近似方法进行处理，感兴趣的读者可参见  Chen et al. (2010); Zymler et al. (2013) 等。

当 $\mathscr{F} \triangleq \mathscr{F}_{WKS}$ 且其中 $K = 1, \underline{p}_1 = \overline{p}_1 = 1$ 时，Hanasusanto et al. (2015) 证明了 (4.20) 等价于如下一系列锥优化约束： 
$$
\begin{array}{lll}
\mathbb{E}ta + \boldsymbol{b}^\top \boldsymbol{\gamma} \ge (1 - \epsilon) \tau,
& \mathbb{E}ta + \boldsymbol{c}_1^\top \boldsymbol{\phi} \le \tau,
& \mathbb{E}ta + \boldsymbol{c}_1^\top \boldsymbol{\psi} \le -y^0(\boldsymbol{x}), \\
\boldsymbol{A}^\top \boldsymbol{\gamma} = \boldsymbol{C}_1^\top \boldsymbol{\phi},
& \boldsymbol{B}^\top \boldsymbol{\gamma} = \boldsymbol{D}_1^\top \boldsymbol{\phi}, \\
\boldsymbol{A}^\top \boldsymbol{\gamma} + \boldsymbol{y}(\boldsymbol{x}) = \boldsymbol{C}_1^\top \boldsymbol{\psi},
& \boldsymbol{B}^\top \boldsymbol{\gamma} = \boldsymbol{D}_1^\top \boldsymbol{\psi}, \\
\mathbb{E}ta \in \mathbb{R}, \boldsymbol{\gamma} \in \mathbb{R}^L, \tau \in \mathbb{R}_+,
&  \boldsymbol{\phi}, \boldsymbol{\psi} \in \mathscr{K}_1^*
\end{array}
$$

其中 $\mathscr{K}_1^*$ 表示 $\mathscr{K}_1$ 的对偶锥 (Ben-Tal and Nemirovski, 2001)。考虑到 $\mathscr{F}_M$ 和 $\mathscr{F}_{DY}$ 均为 $\mathscr{F}_{WKS}$ 的特例，当 $\mathscr{F} \triangleq \mathscr{F}_M$ (其中支撑集为 (4.10)的形式) 或 $\mathscr{F} \triangleq \mathscr{F}_{DY}$ 时，(4.4)可等价转化成锥优化约束形式。

作为 Hanasusanto et al. (2015) 的扩展，Xie and Ahmed (2018) 考虑了更一般的模糊集，研究了分布鲁棒独立机会约束规划和联合机会约束规划的等价凸优化形式。对于考虑均值、散度上界、支撑集的一类模糊集，Hanasusanto et al. (2017) 研究了其分布鲁棒联合机会约束规划的计算复杂度及求解方法。

当 (4.20) 中模糊集取 $\mathscr{F} \triangleq \mathscr{F}_{KL}$ 时，Jiang and Guan (2016) 的研究表明，(4.20) 等价于
$$
\begin{equation}
\label{eq.dro.cc.kl.eq}
\displaystyle \hat{\mathbb{P}}[y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top\tilde{\boldsymbol{\varepsilon}} \le 0] \ge 1 - \bar{\epsilon}.
\end{equation} \tag{4.23}
$$

其中，
$$
\bar{\epsilon} \triangleq 1 - \inf_{t \in (0, 1)} \frac{e^{-\theta}t^{1 - \epsilon} - 1}{t - 1}.
$$
由此可见，它与随机规划中基于采样平均近似 (sample average approximation) 的机会约束式 (4.24) 相比，仅仅是具有不同的概率界 $\bar{\epsilon}$ 而已。此外，对于一般的$\phi$散度形式下分布鲁棒联合机会约束规划的处理方法，可详见  Jiang and Guan (2016)。
$$
\begin{equation}
\label{eq.dro.cc.saa}
\displaystyle \hat{\mathbb{P}}[y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top\tilde{\boldsymbol{\varepsilon}} \le 0] \ge 1 - \epsilon. \tag{4.24}
\end{equation}
$$
作为常用技巧，(4.24)可通过引入0-1辅助决策变量的方法等价转换为
$$
\begin{equation}
\label{eq.dro.cc.saa.eq}
\begin{array}{ll}
\displaystyle y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top\hat{\boldsymbol{\varepsilon}}_\omega\le M_0 (1 - z_\omega), & \forall \omega \in [N], \\
\displaystyle \frac{\sum_{\omega \in [N]} z_\omega}{N} \ge 1 - \epsilon, \\
\boldsymbol{z} \in \{0, 1\}^N.
\end{array}
\end{equation}
$$

其中，$M_0$ 为一个足够大的实数。观察可知，当$y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top\hat{\boldsymbol{\varepsilon}}_\omega\le 0$ 时，$z_\omega$ 可取值 1，代表该约束在第 $\omega$ 个场景中成立，否则不得不取值 0。

当 $\mathscr{F} \triangleq \mathscr{F}_W$ 时，Chen et al. (2018); Xie (2019) 讨论了如何处理分布鲁棒（独立和联合）机会约束规划问题，他们用不同的方法得到了相同的结论。此时，独立机会约束(4.20)等价于如下混合 0-1 锥优化约束
$$
\begin{array}{ll}
\epsilon N t - \boldsymbol{e}^\top \boldsymbol{s} \ge \theta N \lVert \boldsymbol{y}(\boldsymbol{x}) \rVert_*, \\
\boldsymbol{y}(\boldsymbol{x})^\top \hat{\boldsymbol{\varepsilon}}_\omega - y^0(\boldsymbol{x}) + M_0 z_\omega \ge t - s_\omega, & \forall \omega \in [N], \\
M_0(1 - z_\omega) \ge t - s_\omega, & \forall \omega \in [N], \\
t \in \mathbb{R}, \boldsymbol{z} \in \{0, 1\}^N, \boldsymbol{s} \in \mathbb{R}^{N}.
\end{array}
$$
其中 $\boldsymbol{e}$ 代表长度为 $N$、元素全为 1 的向量，$\lVert \cdot \rVert_*$ 为 (4.17)中 $\lVert \cdot \rVert$ 对应的对偶范数 (Boyd et al., 2004)。

## 4.3 分布鲁棒线性优化（Distributionally robust linear optimization）

现实世界中很多优化问题可建模或近似为线性规划问题。线性约束不仅本身可描述许多现实问题的资源约束，而且可建模或近似更复杂的资源约束。作为抛砖引玉，本节讨论如何处理分布鲁棒线性优化约束式 (4.5})，读者可进而自行研究更一般也更难处理的非线性约束，例如：
$$
\begin{equation}
\label{dro.con.pl}
\mathbb{E}_{\mathbb{P}} \big[\max_{k \in [K]} \{y_k^0(\boldsymbol{x}) + \boldsymbol{y}_k(\boldsymbol{x})^\top\tilde{\boldsymbol{\varepsilon}}\}\big] \le 0, \quad \forall \mathbb{P} \in \mathscr{F},
\end{equation} \tag{4.26}
$$
其左端项为关于 $\boldsymbol{x}$ 和 $\tilde{\boldsymbol{\varepsilon}}$ （各自）的分段线性凸函数；此形式出现在许多管理科学问题中，如库存管理 (See and Sim, 2010; Mamani et al., 2017)、预约调度  (Mak et al., 2015; Kong et al., 2013; Qi, 2017)、带时间窗的车辆路径问题  (Zhang et al., 2019)，等等。

接下来探讨如何处理分布鲁棒线性优化约束式 (4.5)。易知，它等价于
$$
\begin{equation}
\label{dro.con.lo.eq}
\sup_{\mathbb{P} \in \mathscr{F}} \mathbb{E}_\mathbb{P}[y^0(\boldsymbol{x}) + \boldsymbol{y}(\boldsymbol{x})^\top \tilde{\boldsymbol{\varepsilon}}] \le 0.
\end{equation} \tag{4.27}
$$
处理该约束的关键是考察左端项中优化问题
$$
\begin{equation}
\label{dro.mod}
Z_P(\boldsymbol{x}) = \sup_{\mathbb{P} \in \mathscr{F}} \mathbb{E}_\mathbb{P}[\boldsymbol{y}(\boldsymbol{x})^\top \tilde{\boldsymbol{\varepsilon}}]
\end{equation} \tag{4.28}
$$

的对偶问题。注意，该优化问题中 $\boldsymbol{x}$ 被视为给定参数，而概率分布 $\mathbb{P}$ 才是决策变量。抽象地，(4.28) 的对偶问题形式为
$$
\begin{equation}
\label{dro.mod.dual}
Z_D(\boldsymbol{x}) = \inf_{\boldsymbol{p} \in \mathscr{P}(\boldsymbol{x})} f(\boldsymbol{p}).
\end{equation}
$$
其中 $\boldsymbol{p}$ 为对偶决策变量，$\mathscr{P}(\boldsymbol{x})$ 为其可行域，$f(\boldsymbol{p})$ 为目标函数，$\boldsymbol{y}(\boldsymbol{x})$ 作为参数被包含于 $\mathscr{P}(\boldsymbol{x})$ 中。在某些条件下，强对偶定理对此成立，则 $Z_P = Z_D$。于是，(4.27) 等价于
$$
\begin{equation}
\label{eq.dro.lo.eq}
\begin{array}{ll}
y^0(\boldsymbol{x}) + f(\boldsymbol{p}) \le 0, \\
\boldsymbol{p} \in \mathscr{P}(\boldsymbol{x}).
\end{array} \tag{4.30}
\end{equation}
$$

因此，技术上主要关注如何在取不同模糊集 $\mathscr{F}$ 的情况下求解 (4.28) 的对偶问题并证明强对偶定理成立。

当 $\mathscr{F} \triangleq \mathscr{F}_{MV}$ 或 $\mathscr{F} \triangleq \mathscr{F}_{MVS}$ 时，由于已知 $\tilde{\boldsymbol{\varepsilon}}$ 的均值 $\boldsymbol{\mu}$，故 (4.28)等价于 $Z_P(\boldsymbol{x}) = \boldsymbol{y}(\boldsymbol{x})^\top \boldsymbol{\mu}$。Popescu (2007) 针对 $\mathscr{F} \triangleq \mathscr{F}_{MV}$ 且 (4.28) 目标函数变为某一类非线性函数的情形，研究了其等价模型与求解方法。

当 (4.28)中 $\mathscr{F} \triangleq \mathscr{F}_{DY}$ 时，则在某些技术性条件下，Delage and Ye (2010)  推导出其对偶问题 (4.29) 的具体形式：
$$
\begin{equation}
\label{dro.mod.dy.dual}
\begin{array}{rll}
\displaystyle Z_D(\boldsymbol{x}) = \min_{\boldsymbol{Q}, \boldsymbol{q}, r, t} & r + t \\
\mbox{s.t.} & r \ge \boldsymbol{y}(\boldsymbol{x})^\top \boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^\top \boldsymbol{Q} \boldsymbol{\varepsilon} - \boldsymbol{\varepsilon}^\top \boldsymbol{q},  \forall \boldsymbol{\varepsilon} \in \Xi, \\
& t \ge (\gamma_2 \boldsymbol{\Sigma} + \boldsymbol{\mu} \boldsymbol{\mu}^\top) \bullet \boldsymbol{Q} + \boldsymbol{\mu}^\top \boldsymbol{q} + \sqrt{\gamma_1} \lVert \boldsymbol{\Sigma}^{1/2} (\boldsymbol{q} + 2 \boldsymbol{Q} \boldsymbol{\mu}) \rVert, \\
& \boldsymbol{Q} \succeq \boldsymbol{0},
\end{array} \tag{4.31}
\end{equation}
$$
其中 ``$\bullet$'' 表示矩阵间的弗罗贝尼乌斯内积。注意到 (4.31) 中的第一个约束实则为（传统）鲁棒优化约束，因此求解 $\mathscr{F}_{DY}$ 模糊集下的分布鲁棒优化问题 (4.28) 等价于求解鲁棒优化问题 (4.31)，而上一章已讲述如何求解鲁棒优化问题。

当 (4.28) 中 $\mathscr{F} \triangleq \mathscr{F}_{WKS}$ 时，则在某些技术性条件下，Wiesemann et al. (2014)推导出其对偶问题 (4.29) 的具体形式：
$$
\begin{equation}
\label{dro.mod.wks.dual}
\begin{array}{rll}
\displaystyle Z_D(\boldsymbol{x}) = \min_{\boldsymbol{\mathbb{E}ta}, \boldsymbol{\eta}, \boldsymbol{\lambda}, \boldsymbol{\phi}} & \displaystyle \boldsymbol{b}^\top \boldsymbol{\mathbb{E}ta} + \sum_{k \in [K]} \overline{\boldsymbol{p}}_k \boldsymbol{\eta}_k - \underline{\boldsymbol{p}}_k \boldsymbol{\lambda}_k \\
\mbox{s.t.} & \displaystyle \boldsymbol{c}_k^\top \boldsymbol{\phi}_k \le \sum_{k' \in \mathcal{A}(k)} (\boldsymbol{\eta}_{k'} - \boldsymbol{\lambda}_{k'}), & \forall k \in [K], \\
& \boldsymbol{C}_k^\top \boldsymbol{\phi}_k + \boldsymbol{A}^\top \boldsymbol{\mathbb{E}ta} = \boldsymbol{y}(\boldsymbol{x}), & \forall k \in [K], \\
& \boldsymbol{D}_k^\top \boldsymbol{\phi}_k + \boldsymbol{B}^\top \boldsymbol{\mathbb{E}ta} = \boldsymbol{0}, & \forall k \in [K], \\
& \boldsymbol{\phi}_k \in \mathscr{K}_k^*, & \forall k \in [K],
\end{array} \tag{4.32}
\end{equation}
$$

其中，$\mathcal{A}(k) \triangleq \{k' \in [K]~|~ \Xi_{k'} \mbox{严格包含于} \Xi_{k} \}$，$\mathscr{K}_k^*$ 表示 $\mathscr{K}$ 的对偶锥 。

当 (4.28)中 $\mathscr{F} \triangleq \mathscr{F}_{KL}$ 时，则在某些技术性条件下，Hu and Hong (2013)推导出其对偶问题 (4.29) 的具体形式：
$$
\begin{equation}
\label{dro.mod.kl.dual}
\begin{array}{rll}
\displaystyle Z_D(\boldsymbol{x}) = \min_{\alpha \ge 0} & \alpha \log \mathbb{E}_{\hat{\mathbb{P}}}[e^{\boldsymbol{y}(\boldsymbol{x})^\top \tilde{\boldsymbol{\varepsilon}}^\dagger / \alpha}] + \alpha \theta.
\end{array} \tag{4.33}
\end{equation}
$$

其中的目标函数为凸函数，因此可用内点法  (Ben-Tal et al., 2013) 或分段线性函数逼近(Long and Qi, 2014)等方法进行处理。

当 (4.28)中 $\mathscr{F} \triangleq \mathscr{F}_{W}$ 且其中的支撑集为 $\Xi \triangleq \{\boldsymbol{\varepsilon} \in \mathbb{R}^I~|~\boldsymbol{C}\boldsymbol{\varepsilon} \le \boldsymbol{d}\}$ 时，则在某些技术性条件下，Esfahani and Kuhn (2018) 推导出其对偶问题 (4.29) 的具体形式：
$$
\begin{equation}
\label{dro.mod.w.dual}
\begin{array}{rll}
\displaystyle Z_D(\boldsymbol{x}) = \inf_{\lambda, \boldsymbol{s}, \boldsymbol{\gamma}} &\displaystyle \lambda \theta + \frac{1}{N} \sum_{\omega \in [N]} s_\omega \\
\mbox{s.t.} & \boldsymbol{y}(\boldsymbol{x})^\top \hat{\boldsymbol{\varepsilon}}_\omega + \boldsymbol{\gamma}_\omega^\top (\boldsymbol{d} - \boldsymbol{C} \hat{\boldsymbol{\varepsilon}}_\omega) \le s_\omega, & \forall \omega \in [N], \\
& \lVert \boldsymbol{C}^\top \boldsymbol{\gamma}_\omega - \boldsymbol{y}(\boldsymbol{x}) \rVert_* \le \lambda, & \forall \omega \in [N], \\
& \boldsymbol{\gamma}_\omega \ge \boldsymbol{0} & \forall \omega \in [N].
\end{array} \tag{4.34}
\end{equation}
$$


将以上各结果嵌入到 (4.30)即可得到 (4.27)的等价形式，进而求解鲁棒线性规划问题。 

事实上，Popescu (2007); Delage and Ye (2010); Wiesemann et al. (2014); Esfahani and Kuhn (2018) 的研究均解决了 (4.26) 的等价转化问题，而上述结论针对 (4.27)，只是(4.26) 中 $K = 1$ 时的特例，感兴趣的读者可细读他们的论文。

## 参考文献

**Ben-Tal, Aharon and Arkadi Nemirovski,** Lectures on modern convex optimization: analysis, algorithms, and engineering applications, Vol. 2, Siam, 2001.

**__ , Dick Den Hertog, Anja De Waegenaere, Bertrand Melenberg, and Gijs Rennen,** “Robust solutions of optimization problems affected by uncertain probabilities,” Management Science, 2013, 59 (2), 341–357. 

**Bertsimas, Dimitris and Ioana Popescu,** “Optimal inequalities in probability theory: A convex optimization approach,” SIAM Journal on Optimization, 2005, 15 (3), 780–804.

**__ , Melvyn Sim, and Meilin Zhang**, “Adaptive distributionally robust optimization,” Management Science, 2019, 65 (2), 604–618. 

**Boyd, Stephen, Stephen P Boyd, and Lieven Vandenberghe**, Convex optimization, Cambridge university press, 2004. 

**Charnes, Abraham and William W Cooper,** “Chance-constrained programming,” Management Science, 1959, 6 (1), 73–79. 

**Chen, Wenqing and Melvyn Sim,** “Goal-driven optimization,” Operations Research, 2009, 57 (2), 342–357. 

**__ , __ ,** **Jie Sun, and Chung-Piaw Teo**, “From CVaR to uncertainty set: Implications in joint chance-constrained optimization,” Operations Research, 2010, 58 (2), 470–485.

**Chen, Zhi, Daniel Kuhn, and Wolfram Wiesemann,** “Data-driven chance constrained programs over Wasser- stein balls,” Available at arXiv:1809.00210, 2018. 

**Delage, Erick and Yinyu Ye,** “Distributionally robust optimization under moment uncertainty with application to data-driven problems,” Operations Research, 2010, 58 (3), 595–612. 

**Esfahani, Peyman Mohajerin and Daniel Kuhn,** “Data-driven distributionally robust optimization using the Wasserstein metric: Performance guarantees and tractable reformulations,” Mathematical Programming, 2018, 171 (1-2), 115–166. 

**Gao, Rui and Anton J Kleywegt,** “Distributionally robust stochastic optimization with Wasserstein distance,” Available at arXiv:1604.02199, 2016. 

**Ghaoui, Laurent El, Maksim Oks, and Francois Oustry,** “Worst-case value-at-risk and robust portfolio optimization: A conic programming approach,” Operations Research, 2003, 51 (4), 543–556. 

**Hanasusanto, Grani A, Vladimir Roitch, Daniel Kuhn, and Wolfram Wiesemann,** “A distributionally robust perspective on uncertainty quantification and chance constrained programming,” Mathematical Programming, 2015, 151 (1), 35–62. 

**__ , __ , __, and __ ,** “Ambiguous joint chance constraints under mean and dispersion information,” Operations Research, 2017, 65 (3), 751–767. 

**Hu, Zhaolin and Lj Hong**, “Kullback-Leibler Divergence -Constrained Distributionally Robust Optimization,” Optimization Online, 2013, (2), 1–34. 

**Jiang, Ruiwei and Yongpei Guan,** “Data-driven chance constrained stochastic program,” Mathematical Pro- gramming, 2016, 158 (1-2), 291–327. 

**Kong, Qingxia, Chung-Yee Lee, Chung-Piaw Teo, and Zhichao Zheng,** “Scheduling arrivals to a stochastic service delivery system using copositive cones,” Operations Research, 2013, 61 (3), 711–726. 

**Long, Daniel Zhuoyu and Jin Qi,** “Distributionally robust discrete optimization with Entropic Value-at-Risk,” Operations Research Letters, 2014, 42 (8), 532–538. 

**Mak, Ho-Yin, Ying Rong, and Jiawei Zhang,** “Appointment scheduling with limited distributional informa- tion,” Management Science, 2015, 61 (2), 316–334. 

**Mamani, Hamed, Shima Nassiri, and Michael R Wagner,** “Closed-form solutions for robust inventory management,” Management Science, 2017, 63 (5), 1625–1643. 

**Natarajan, Karthik, Chung Piaw Teo, and Zhichao Zheng,** “Mixed 0-1 linear programs under objective uncertainty: A completely positive representation,” Operations Research, 2011, 59 (3), 713–728. 

**__ , Dessislava Pachamanova, and Melvyn Sim,** “Constructing Risk Measures from Uncertainty Sets,” Oper- ations Research, 2009, 57 (5), 1129–1141. 

**Nemirovski, Arkadi and Alexander Shapiro,** “Convex approximations of chance constrained programs,” SIAM Journal on Optimization, 2006, 17 (4), 969–996. 

**Popescu, Ioana,** “Robust mean-covariance solutions for stochastic optimization,” Operations Research, 2007, 55 (1), 98–112. 

**Qi, Jin,** “Mitigating delays and unfairness in appointment systems,” Management Science, 2017, 63 (2), 566–583. 

**See, Chuen-Teck and Melvyn Sim,** “Robust approximation to multiperiod inventory management,” Operations Research, 2010, 58 (3), 583–594. 

**Wiesemann, Wolfram, Daniel Kuhn, and Melvyn Sim,** “Distributionally robust convex optimization,” Oper- ations Research, 2014, 62 (6), 1358–1376. 

**Xie, Weijun**, “On distributionally robust chance constrained programs with Wasserstein distance,” Mathematical Programming, 2019, pp. 1–41. 

**__ and Shabbir Ahmed,** “On deterministic reformulations of distributionally robust joint chance constrained optimization problems,” SIAM Journal on Optimization, 2018, 28 (2), 1151–1182. 

**Zhang, Yu, Roberto Baldacci, Melvyn Sim, and Jiafu Tang,** “Routing optimization with time windows under uncertainty,” Mathematical Programming, 2019, 175 (1-2), 263–305.

**Zhao, Chaoyue and Yongpei Guan,** “Data-driven risk-averse stochastic optimization with Wasserstein metric,” Operations Research Letters, 2018, 46 (2), 262–267. 

**Zymler, Steve, Daniel Kuhn, and Berç Rustem,** “Distributionally robust joint chance constraints with second- order moment information,” Mathematical Programming, 2013, 137 (1-2), 167–198. 
