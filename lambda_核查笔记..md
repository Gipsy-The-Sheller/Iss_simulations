# λ 统计量：原文口径核查与延伸核验

> 核查对象：Xu, Z.-J. (2026) *What the index of substitution saturation measures, and why its critical value is a power threshold*，§2.15、Proposition 23、Table 3、§4.5、摘要。
> 对照原始出处：Jermiin & Misof (2020), *Measuring Historical and Compositional Signals in Phylogenetic Data*, bioRxiv 10.1101/2020.01.03.894097。
> 核验脚本：`lhcore.py`（谱分解 + Felsenstein pruning）、`lam.py`（本节全部数值）。
>
> **体例说明**：「原文」= 预印本文字，标注章节；「核验」= 本文复算或复现；「推断」= 我的延伸，非原文主张。三者严格分开。

---

## 0. 一句话结论

原文的批评**只对 $\bar\lambda$ 成立，对逐对 $\lambda_{ij}$ 不成立**。而且原文自身的 Lemma 2 已经蕴含了 λ 的一个正面性质（δ_ran 就是任意可逆模型的饱和上界），但它没有把这层说出来。

---

## 1. 原文口径：逐字核对

### 1.1 定义（§2.15，式 20）

> λ is defined per pair as
> $$\lambda_{ij}=\frac{\delta_{obs}(ij)}{\delta_{ran}(ij)},\qquad \delta_{ran}(ij)=1-\sum_k f_{ik}f_{jk},$$
> where $f_{ik}$ is the frequency of state k in sequence i. The denominator is the divergence expected of two sequences that have evolved independently to equilibrium.

**要点**：定义是**逐对**的（per pair），且 $f_{ik}$、$f_{jk}$ **分开取**——两条序列各自有自己的组成。这一点后面是关键。

### 1.2 Proposition 23（§2.15）

> **Proposition 23 (λ is the p-distance rescaled).** If every sequence has composition q, then $\lambda_{ij}=p_{ij}/(1-\sum_k q_k^2)$ for every pair, and therefore $\bar\lambda=\bar p/(1-\sum_k q_k^2)$.

> **Proof.** δ_obs(ij) is by construction the proportion of differing sites between i and j, i.e. $p_{ij}$; δ_ran(ij) reduces to $1-\sum_k q_k^2$ under the stated condition. Average over pairs.

**注意命题的假设**：「If every sequence has composition q」——这是**共享组成**假设。原文自己在上一句说明，这是"the second-order approximation the rest of this paper works in"（全文其余部分采用的二阶近似）。

### 1.3 原文的数值核验（§2.15）

> Table 3 evaluates it on seven alignments spanning $I_{ss}$ from 0.038 to 0.749, and the right-hand side of Proposition 23 reproduces the observed $\bar\lambda$ **to four decimal places** on every one of them. Across the seven, λ correlates with $\bar p$ at r = 0.999 and with $I_{ss}$ at r = 0.998.

### 1.4 原文的定性判断（§2.15）

> Whatever else λ may be — and it has real advantages, **being defined pairwise and permitting non-stationary conditions** — it is not an independent measurement of saturation: it is the uncorrected proportion of differing sites in different units, and it inherits everything Theorem 13 establishes about that quantity.

**这是最需要精确读的一句。** 原文明确承认了 λ 的两个优点：**逐对定义**、**允许非平稳**。然后才说它"不是独立测量"。

### 1.5 摘要的口径

> The same holds for λ, the principal tree-free alternative in use: under a shared base composition it is exactly the mean p-distance rescaled by $1/(1-\sum_k q_k^2)$, which reproduces the computed values to four decimal places across seven alignments spanning the index from 0.038 to 0.749.

摘要同样带限定语"under a shared base composition"。

### 1.6 方法部分（§4.5）

> λ of Jermiin & Misof [2020], defined pairwise by Eq. (20). Both are computed here on the same alignments as $I_{ss}$, from their published definitions.

### 1.7 表 3 的呈现方式

> Table 3: ... **λ is the mean over all pairs.** It tracks $I_{ss}$ closely, as Proposition 23 requires, and the last column is its ratio to the mean p-distance, which is $1/(1-\sum_k q_k^2)$ and therefore nearly constant.

---

## 2. 核查一：「λ 是 JC69 特有的吗？」——不是

**原文未主张 λ 与 JC69 有任何关系。** 通篇 λ 只在 §2.15、摘要、§4.5 出现，均无任何模型假设。Proposition 23 的推导只用到了 δ_obs = p 和共享组成，**不涉及 Q**。

### 核验：δ_ran 就是任意可逆模型的饱和上界

由本文 Lemma 2：

$$\wp=\Big(1-\sum_j\pi_j^2\Big)-\sum_{\mu\neq0}w_\mu\xi_\mu,\qquad \xi_\mu=\langle e^{\mu d}\rangle$$

令 $\xi_\mu\to0$（饱和极限），得

$$\wp(\infty)=1-\sum_j\pi_j^2$$

而共享组成 $q$ 下 $\delta_{ran}=1-\sum_k q_k^2$。**两者是同一个量。**

> **推断**：原文 Lemma 2 已经蕴含「δ_ran = 任意可逆模型下 p-distance 的饱和上界」，但 §2.15 没有把这层指出来。这是原文的一处**未展开**（不是错误）。

### 核验：0.75 只是一个点

$1-\sum_k q_k^2$ 随字母表与组成变化：

| 组成 | K | 上界 | 与 0.75 之差 |
|---|---|---|---|
| 核苷酸均匀（JC69 隐含） | 4 | 0.7500 | — |
| (.40,.30,.20,.10) | 4 | 0.7000 | −0.050 |
| (.40,.40,.10,.10) | 4 | 0.6600 | −0.090 |
| (.45,.45,.05,.05) | 4 | 0.5900 | −0.160 |
| 氨基酸均匀 | 20 | 0.9500 | **+0.200** |
| 氨基酸实际频度 | 20 | 0.9315 | **+0.182** |

**氨基酸那两行是最危险的**：20 态下上界约 0.93–0.95，用 0.75 当阈值会把远未饱和的蛋白比对**反向误判**为饱和。

### 核验：非 JC、非均匀模型下 λ<1 与内点 MLE 的等价性

GTR（随机交换率）+ π=(.40,.30,.20,.10)，谱 $\lambda=(-0.8451,-1.5259,-1.7953)$，完全非 JC69、非均匀：

| $t_{true}$ | $\hat p$ | 上界 | λ | 状态 | $\hat t$ | $|\ell'(\hat t)|$ |
|---|---|---|---|---|---|---|
| 0.05 | 0.0472 | 0.700 | 0.067 | interior | 0.050 | 3e−12 |
| 0.30 | 0.2170 | 0.700 | 0.310 | interior | 0.310 | 7e−13 |
| 1.00 | 0.4125 | 0.700 | 0.589 | interior | 1.008 | 4e−13 |
| 2.00 | 0.5192 | 0.700 | 0.742 | interior | 2.012 | 2e−13 |
| 3.00 | 0.5736 | 0.700 | 0.819 | interior | 3.078 | 3e−14 |
| 5.00 | 0.6268 | 0.700 | 0.895 | interior | 5.212 | 1e−14 |
| 8.00 | 0.6621 | 0.700 | 0.946 | interior | 8.310 | 1e−14 |

> **推断**：λ<1 ⟺ 存在内点 MLE，成立条件是**可逆 + 平稳**，不是 JC69。

### 核验：从原文自己的 Table 3 反推

用表 3 的 $\lambda/\bar p$ 列反解 $1-\sum q_k^2$ 与 $\bar p$：

| 数据集 | $I_{ss}$ | λ | λ/$\bar p$ | 上界 | $\bar p$ | **$\bar p/0.75$** | **λ** |
|---|---|---|---|---|---|---|---|
| COI 第二位 | 0.038 | 0.041 | 1.404 | 0.712 | 0.029 | 0.039 | 0.041 |
| COI 第一+二位 | 0.087 | 0.096 | 1.344 | 0.744 | 0.071 | 0.095 | 0.096 |
| 鸟类完整列 | 0.119 | 0.108 | 1.345 | 0.744 | 0.080 | 0.107 | 0.108 |
| COI 第一位 | 0.137 | 0.151 | 1.339 | 0.747 | 0.113 | 0.150 | 0.151 |
| 鸟类信息位点 | 0.268 | 0.253 | 1.338 | 0.747 | 0.189 | 0.252 | 0.253 |
| COI 全部位点 | 0.270 | 0.294 | 1.350 | 0.741 | 0.218 | 0.290 | 0.294 |
| **COI 第三位** | 0.715 | 0.749 | **1.462** | **0.684** | 0.512 | **0.683** | **0.749** |

**第三位密码子这一行是全文最有说服力的例子**：它的组成偏斜使上界降到 **0.684**，于是

- 朴素规则 $\bar p/0.75=0.683$ → 「还差 32% 才饱和」
- λ = 0.749 → 「历史信号已损失 75%」

两者相差 0.066 绝对值（相对差约 10%）。**这不是我的构造，是原文表 3 自带的数字。**

---

## 3. 核查二：λ 相对 $\bar p$ 的增量在哪

### 3.1 在共享组成下，对 $\bar\lambda$ 严格为零

Proposition 23 说的是常数重标度：$\bar\lambda=\bar p/(1-\sum q_k^2)$，而 $1-\sum q_k^2$ 只依赖 pooled 组成，**不依赖任何树或分歧信息**。所以：

> 作为**单一标量**，$\bar\lambda$ 与 $\bar p$ 携带的信息严格相同。原文此处的批评成立。

### 3.2 逐对结构 —— 真正的增量

原文在 §2.15 承认 λ「being defined pairwise」，但表 3 只报了 $\bar\lambda$（"λ is the mean over all pairs"），**把这一层压掉了**。

**核验**（六类群，一条长枝 Z–W=1.6 且 F 挂在其末端，其余短枝 0.05，Γ(α=1)，4000 位点）：

| 全局量 | 值 |
|---|---|
| $I_{ss}$ | 0.259 |
| $\bar p$ | 0.266 |
| $\bar\lambda$ | 0.354 |

逐对 $\lambda_{ij}$：

| 含 F 的五对 | 其余十对 |
|---|---|
| B–F 0.725, D–F 0.720, A–F 0.715, C–F 0.711, E–F 0.709 | 0.114 – 0.218 |

**$\bar\lambda=0.354$ 是一个"平均"出来的中间值，看不出任何异常；逐对矩阵直接把 F 那一条枝点名。** 这正是 RedundancyHeatMapper（Misof, 2019）设计成画 λ 矩阵热图的原因。

> **推断**：原文把 λ 压成均值再比较，恰好抹掉了它唯一真正的信息。批评「λ 是换单位的 p-distance」对 $\bar\lambda$ 成立，对 $\lambda_{ij}$ 不成立。

### 3.3 非平稳定义

$f_{ik}$、$f_{jk}$ 分开取，两条序列组成不同时 δ_ran 仍有定义。

> **推断（重要限定）**：此时 $\delta_{ran}=1-\sum_k f_{ik}f_{jk}$ 不再是模型导出的饱和上界，而是**经验随机化基线**，§2 的等价性（λ<1 ⟺ 内点 MLE）**不再成立**。这是 λ 的"非平稳优势"的代价，原文未讨论。

---

## 4. 核查三：λ 的一个硬伤（原文未提，Jermiin & Misof 原文亦未处理）

λ 用的是**未校正的** $\hat p$，没有标准误。

**核验**：真值恰好落在上界时，$\Pr(\hat p\ge\text{上界})$：

| L | 200 | 500 | 1000 | 5000 | 20000 |
|---|---|---|---|---|---|
| 假阳率 | 0.473 | 0.483 | 0.488 | 0.495 | 0.497 |

**序列长度再大也救不了** —— 恒在 50% 附近。

反过来看功效（真 $\lambda=0.99$ 时误判为饱和的概率）：

| 真 λ | L=200 | L=1000 | L=5000 | L=20000 |
|---|---|---|---|---|
| 0.90 | 0.016 | 0.000 | 0.000 | 0.000 |
| 0.95 | 0.130 | 0.008 | 0.000 | 0.000 |
| 0.98 | 0.310 | 0.162 | 0.016 | 0.000 |
| 0.99 | 0.389 | 0.305 | 0.138 | 0.015 |

> **推断**：λ<1 是一个**没有标准误的检验**。它缺的正是 $(\hat t,\ \ell'')$ 方案给的东西——曲率。$\hat t$ 的 SE 直接回答"距上界还差几个标准误"，λ 只给一个裸比值。

---

## 5. 与作者私下表述的对接

作者的判据：算目标枝的剖面似然 $\ell(t_b)$，看 $\ell'$ 是否**过零**（内点解）而非"山峰落在模型范围外"。

| 判据 | 形式 | 等价表述 |
|---|---|---|
| 一阶导过零 | $\exists t^*:\ell'(t^*)=0,\ \ell''<0$ | — |
| λ 阈值 | $\lambda<1$ | $\hat p<1-\sum q_k^2$ |
| SatuTe 的 z 分数 | 相干系数显著性 | 两子树间分歧是否可辨 |

**核验**：两序列 JC69 下 $\ell'(t)\propto 4n(\tfrac43\hat p-1)$，故内点 MLE ⟺ $\hat p<0.75$ ⟺ λ<1。多分类时"表观 p-distance"换成"两子树间分歧"，就是 SatuTe 在测的东西。

> **推断**：这三种表述是同一个"实际可识别性"问题的三种形式。第一种最优，因为同时给出「能不能估」（过零）和「估得多差」（$\ell''$ → SE）；λ<1 是内点解的**必要不充分**条件（在两序列 JC 下恰好等价），且不给曲率。

---

## 6. 需要修正的前期表述

| 先前说法 | 修正 |
|---|---|
| 「λ 是 JC69-specific？」 | **否**。δ_ran 模型无关；0.75 才是均匀组成那一点。 |
| 「原文说 λ 只是换单位的 p-distance」 | 原文限定在**共享组成**下、且对象是 $\bar\lambda$。它同时明确承认 λ「defined pairwise and permitting non-stationary conditions」。 |
| 「λ 相对 p̄ 的增量为零」 | 对 $\bar\lambda$ 成立；对逐对 $\lambda_{ij}$ **不成立**（见 §3.2 核验）。 |

---

## 7. 实用建议

1. **不要用 0.75**。用 $1-\sum_k q_k^2$，按你自己的 pooled 组成算。氨基酸数据尤其重要（上界 0.93–0.95）。
2. **看逐对矩阵，不看均值**。用 SatuRation v1.0 + RedundancyHeatMapper，或自己算。
3. **λ 只当筛子用，不当检验用**。它没有 SE，在上界附近功效约 50%。
4. **要定量就用似然**。$(\hat t,\ell'')$ 三元组同时给出可估性与精度，严格强于 λ。
5. **别把 λ 和 $I_{ss}$ 混用定位**：$I_{ss}$ 是全局标量无法定位；λ 的矩阵可以。

---

## 附：复现

```bash
python3 lam.py      # 上界表、GTR 等价性、噪声功效、逐对定位
python3 lamplot.py  # 图：上界随偏斜度变化 / λ 检验功效
```

依赖：`numpy`, `scipy`, `matplotlib`；中文字体 WenQuanYi Micro Hei。
