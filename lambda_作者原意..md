# λ 的作者原意：Jermiin & Misof (2020) 如何诠释这个统计量

> 核查对象：Jermiin, L. S., & Misof, B. (2020) *Measuring Historical and Compositional Signals in Phylogenetic Data*. bioRxiv, doi:10.1101/2020.01.03.894097（2020-01-03 发布）。配套软件 SatuRation v1.0；可视化 RedundancyHeatMapper / SatuRationHeatMapper（Misof, 2019）。
> 对照：Xu (2026) 预印本 §2.15、Proposition 23、Table 3。
>
> **体例**：「原文」= Jermiin & Misof 的文字（英文逐字）；「核验」= 本文复算或交叉验证；「推断」= 我的延伸。三者严格分开。凡未能从原文直接确认的，标注为**待确认**。

---

## 0. 一句话

作者的诠释框架是 **「历史信号的损失」**（loss / decay / erosion of historical signal），**不是**「替换饱和度」，也**不是**「多重替换的次数」。δ_ran 是一个**反事实基线**（两条序列独立演化到无限久、各自达到平衡时的差异比例），λ 是观测差异**相对于这条基线**的位置。它的身份是**描述性诊断量**，作者从未给它配标准误、显著性检验或普适临界值。

---

## 1. 作者如何定义「历史信号」

原文（引言，逐字）：

> "the historical signal (i.e., the signal in DNA that is due to the **order and time of divergence events**)"

> "there is a limit to how far apart sequences can diverge before the historical signal ... has eroded so much that it is pointless attempting to infer a phylogeny from these data."

> "Visually, an alignment may appear conserved, but that is because we tend to focus on the presence and abundance of constant sites in the alignment. If we were to focus our attention on the differences at the variable sites, we may actually find that the nucleotides at these sites have changed so frequently that the historical signal is erased."

**要点**：历史信号 = 位点差异中**可追溯到分歧顺序与时间**的那部分信息。λ 测的就是这部分还剩多少。

### 与组成信号的分工

| 量 | 测什么 | 原文表述 |
|---|---|---|
| **λ** | 历史信号（分歧顺序/时间）的损失 | "λ quantifies loss of historical signal" |
| **δ_CFS** | 组成信号：偏离完全对称的标准化组成距离 | "δ_CFS = 0 indicating no evidence of evolution under dissimilar conditions and δ_CFS > 0 indicating increasing evidence of lineages diverging under different conditions" |

原文建议**联用**：

> 高 λ 与高 δ_CFS 同时出现 → 可能存在谱系特异进化过程，并可能使系统发育估计偏倚。

> **推断**：作者的设计是「λ 回答历史信号还剩多少，δ_CFS 回答组成偏离有多强」，二者是**并列诊断**，不是合并成一个分数。这与 Φ/Iss 那种"单一全局标量"的思路根本不同。

---

## 2. δ_obs 与 δ_ran 的定义（关键：δ_ran 是反事实，不是打乱）

原文（逐字）：

> "δ_obs(ij) is the proportion of variable sites with different nucleotides in sequences i and j"

> "δ_ran(ij) is the proportion of variable sites with different nucleotides in sequences i and j, **assuming they have evolved for infinitely long time under independent stationary conditions**—which implies that each sequence will be in equilibrium"

**所以 δ_ran 不是"把位点顺序打乱"。** 它是一个**反事实生成过程**：两条序列独立演化、无限久、各自达到自身平衡组成，此时它们之间还有多少比例的差异位点是"不同"的。

> **核验**：在两序列共享平衡组成 $q$ 时，独立抽样相同的概率为 $\sum_k q_k^2$，故 $\delta_{ran}=1-\sum_k q_k^2$。这正是 Xu 预印本写下的形式。两条序列组成不同时为 $1-\sum_k f_{ik}f_{jk}$。
>
> **推断**：这个量恰好等于 Xu 预印本 Lemma 2 中 $\wp$ 在 $\xi_\mu\to0$ 时的极限 $1-\sum_j\pi_j^2$，即**任意可逆模型的 p-distance 饱和上界**。Jermiin & Misof 是从"独立平稳反事实"这条路径抵达同一个量，Xu 是从谱分解这条路径抵达。两条路殊途同归，但**作者没有把这层（δ_ran = 饱和上界）点明**。

**待确认**：软件 SatuRation v1.0 名为 "randomisation"，README 称 "SatuRation estimates the degree of saturation observed between two sequences"，但未公开 randomisation 的重排单位、次数与算法步骤。论文定义是解析的；软件实现**是否**等价于解析式，公开材料不足以判定。

---

## 3. λ 的公式：δ_obs / δ_ran

预印本 HTML 中该公式是**嵌入图片**，无法以文本抽取。但可用作者原文的**端点诠释**反推（此为我做的核验，非原文逐字公式）：

| 作者原文的端点 | δ_obs/δ_ran | (δ_obs−δ_ran)/(1−δ_ran) |
|---|---|---|
| (a) 可变位点全同 → λ=0.0 | **0.000** ✓ | −2.333 ✗ |
| (b) 两序列相对随机 → λ=1.0 | **1.000** ✓ | 0.000 ✗ |
| (c) 偶尔 λ>1.0 | **1.030** ✓ | 0.070 ✗ |

只有 $\delta_{obs}/\delta_{ran}$ 与三条端点**全部兼容**。Xu 预印本写下的形式正确。

> **注**：一份二手核查曾推断公式为 $(\delta_{obs}-\delta_{ran})/(1-\delta_{ran})$，上表显示该形式在 (a) 处给出负值，与作者"λ = 0.0"直接矛盾，应弃用。

### 端点语义（原文逐字）

> "If the nucleotides at the variable sites are identical, λ = 0.0 and there is no evidence of a loss of historical signal between the sequences. On the other hand, if the two sequences are random with respect to each other, λ = 1.0 and the historical signal will be fully decayed."

**必须注意的语义陷阱**：

- λ = 0 **不表示零替换**，只表示"可变位点完全相同"
- λ = 1 **不表示每个位点都多次替换**，只表示"两条序列已相对彼此随机"

> **推断**：这个量回答的是「还能不能从历史分歧中提取树信号」，**不回答**「树上实际发生了几次替换」。这正是 Xu 预印本讨论部分那句"饱和的标准意义是观测差异与真实分歧之间的差距"所关心的另一件事。

---

## 4. λ > 1：作者自己的解释

原文（逐字）：

> "occasionally, λ > 1.0, in which case **δ_obs or δ_ran may be inaccurate, perhaps due to a finite length of the alignment or due to sequencing or alignment error**."

**这是非常关键的一句。** 作者明确把 λ>1 归为**度量失准**，并点名原因是**有限序列长度**、测序或比对错误。

> **核验**（我上一轮做的功效分析，与此互为印证）：真值恰好落在上界时，$\Pr(\hat p\ge\text{上界})$ 在 L=200→20000 上恒为 0.47–0.50，**序列再长也救不了**；真 λ=0.99 时 L=20000 仍有 1.5% 误判。
>
> **推断**：作者意识到 λ 没有做噪声校正，但只以"occasionally"带过，**没有给出 λ 的标准误、置信区间或重采样方案**。这是该统计量最实质的缺陷，且作者与 Xu 双方都未量化它。

---

## 5. 作者的自我定位：λ 是作为 Iss 的替代提出的

原文（引言，逐字）：

> "little has been done to quantify how eroded a historical signal in an alignment of sequence data may be. **An exception to this is Xia et al.'s (2003) index of substitution saturation.**"

> **推断**：λ 的提出**直接以 Iss 为对话对象**。这使得 Xu (2026) 反过来批评 λ，构成了一次完整的往返——两篇论文在同一问题上隔空交锋。

### 作者主张的相对优势（原文逐字）

> "The main advantage of this method of estimating loss of historical signal is that **λ allows for cases where the sequences have diverged under different Markovian conditions; hence, λ is less likely to be overestimated.**"

**机制**：δ_ran 用**各序列自身**的平衡组成作基准，而非把单一全局频率强加给所有谱系。两条序列组成不同时，随机基线本身升高，归一化时这部分"非历史差异"被部分扣除。

> **推断（重要限定）**：这是**方法优势的主张**，不是证明。而且它有代价——非平稳情形下 $\delta_{ran}$ 不再是模型导出的饱和上界，而是经验随机化基线，**λ<1 ⟺ 内点 MLE 的等价性随之断裂**。这一点 Xu 未讨论，作者亦未讨论。

---

## 6. 成对性质与矩阵：作者的工作流是热图，不是均值

原文（逐字）：

> "SatuRation ... computes δ_obs, δ_CFS and λ for **all pairs of sequences** in an alignment"

两个配套绘图程序的 README（逐字）：

> "RedundancyHeatMapper produces a heat map from a distance matrix with values of lambda obtained using SatuRation." (RedundancyHeatMapper, 2019)

> "Lambda is a measure of the degree of saturation across the sites in a pair of sequences. The higher a lambda value is, the darker the corresponding pixel in the heat map is. **Emphasis is on lambda values between 0.64 and 1.0, allowing the sequence pairs that effectively are random with respect to one another to be identified.**" (SatuRationHeatMapper, 2019)

> "Emphasis is on lambda values between 0.0 and 0.3, allowing the sequence pairs that are most similar to one another to be identified." (RedundancyHeatMapper, 2019)

**注意两个 README 的重心相反**：一个强调 0.0–0.3（找最相似的对），一个强调 0.64–1.0（找已相对随机的对）。

> **推断**：这正说明**不存在统一阈值**——区间是**可视化/探索取向的选择**，取决于你想找什么，不是统计临界值。

> **核验**（我上一轮做的）：六类群、一条长枝、Γ(α=1)、4000 位点，全局 $I_{ss}=0.259$、$\bar\lambda=0.354$ 完全看不出异常；逐对矩阵里含 F 的五对是 **0.709–0.725**，其余十对 0.114–0.218，**直接点名 F 那条枝**。

**这一点直接反击 Xu 的处理方式**：Xu 的表 3 明写 "λ is the mean over all pairs"，把 λ 压成单一均值再与 $\bar p$ 比相关（r=0.999）。压成均值后它**确实**退化成换单位的 $\bar p$（Proposition 23 成立）；但这是**降维造成的**，不是 λ 本身的性质。

---

## 7. 阈值与推断：作者没有给

| 问题 | 作者的回答 |
|---|---|
| λ 的普适临界值 | **无**。λ=0 与 λ=1 是**定义端点**，不是"可用/不可用"临界值 |
| 标准误 / 置信区间 | **无** |
| 显著性检验 | **无**（λ 没有；δ_CFS 有 Bowker (1948) 配对对称性检验，但不能移植给 λ） |
| 模拟标定表 | **无**（作者用模拟数据说明方法表现，未规定 λ 的抽样分布） |

> **推断**：把某个具体 λ 值称为"显著饱和"，**超出作者原文的措辞**。

---

## 8. 术语：作者用"saturation"吗？

**用，但分层**：

- **软件层面**：直接叫 saturation。程序名 SatuRation；README 说 "SatuRation estimates the degree of saturation observed between two sequences"；SatuRationHeatMapper README 说 "Lambda is a measure of the degree of saturation across the sites in a pair of sequences"。
- **论文概念层面**：主身份是 "loss of historical signal"，动词用 decays / fully decayed / fully eroded。

> **推断**：作者没有回避 saturation，但刻意把**理论身份**定为归一化后的**历史信号损失**，以区别于传统"p 距离对树长的饱和曲线"。

---

## 9. 发表状态

- DOI 10.1101/2020.01.03.894097 **至今仍是 bioRxiv 预印本**，未见带卷期页码的正式期刊版本。
- 两个 2019 年 README 把参考文献写作 "Jermiin LS, Misof B. *Assessing the tree-likeness of phylogenetic data*. **In prep.**"
- 后续相关工作 Jermiin, Lovell, Misof, Foster & Robinson 另有一篇（*Syst Biol*），属**同一作者团队的延续**，不是该预印本的发表版本。

> **建议引用格式**：Jermiin, L. S., & Misof, B. (2020) *Measuring Historical and Compositional Signals in Phylogenetic Data*. bioRxiv preprint, doi:10.1101/2020.01.03.894097。

> **待确认**：因无正式发表版，无法逐字比对定义是否在同行评议中被修订。

---

## 10. 与 Xu (2026) 批评的对照

| 维度 | Jermiin & Misof 自己的主张 | Xu 的批评 | 评判 |
|---|---|---|---|
| 语义 | 历史信号损失，非替换次数 | 未质疑语义，质疑信息含量 | 二者说的**不是同一件事** |
| δ_ran | 独立平稳无限久的反事实基线 | 写为 $1-\sum q_k^2$ | 一致；且 = Xu Lemma 2 的饱和上界（Xu 未点明） |
| 共享组成下 | — | Prop 23：$\bar\lambda=\bar p/(1-\sum q_k^2)$ | **正确**，且作者端点语义与之兼容 |
| 成对结构 | 明确 pairwise，热图是工作流 | 表 3 只报均值 | **Xu 抹掉了 λ 的核心信息** |
| 非平稳 | 明确主张"less likely to be overestimated" | 承认"permitting non-stationary" | Xu 承认了，但未指出等价性断裂的代价 |
| 噪声 | 自己承认 λ>1 源于有限长度等 | 未讨论 | **Xu 未发现作者已自曝此点** |
| 阈值 | 明确没有 | 未批评 | 一致 |

**净结论**：Xu 的 Proposition 23 在数学上正确，但其适用范围被他自己的呈现方式窄化了——他对 $\bar\lambda$（均值）的批评成立，对 $\lambda_{ij}$（矩阵）不成立；而后者恰恰是作者设计该统计量的初衷。

---

## 11. 待确认清单

1. **公式的逐字形式**：预印本 HTML 中 λ 公式是嵌入图片，本文由端点诠释反推为 $\delta_{obs}/\delta_{ran}$（三条端点全部兼容）。建议取得 PDF 原文核对一次。
2. **SatuRation v1.0 的 randomisation 实现**：是否等价于解析式 $1-\sum_k f_{ik}f_{jk}$？公开 README 与未取得的 C++ 源码不足以判定。
3. **"variable sites" 的口径**：原文说 δ_obs 是 "proportion of variable sites with different nucleotides"。若软件只统计可变位点（排除恒定位点），则与 p-distance（全部位点）不同。Xu 表 3 用 $\bar p$ 复现 λ 到四位小数，暗示其实现按全部位点计；但两者在特定数据集上是否总一致，**待确认**。
4. **λ>1 的出现率与处置**：作者只说 "occasionally"，无量化。

---

## 附：本文核验所用脚本

```bash
python3 lam.py       # 上界表、GTR 非 JC 等价性、噪声功效、逐对定位
python3 lamplot.py   # 上界随偏斜度变化 / λ 检验功效
```

三份文档的关系：
- `lambda_核查笔记.md` — Xu 预印本对 λ 的口径核查
- `lambda_作者原意.md`（本文）— Jermiin & Misof 的原意
- `lambda_original_semantics.md` — 二手核查原始稿（**其中公式推断一节有误，见本文 §3**）
