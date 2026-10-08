# 模拟实验计划（Iss / I-distance 批判项目）

> 建立于 2026-10-08。本文件是模拟部分的施工图；数学部分的讨论在 `Iss_JC69_star_proof (1).md`
> 与 `Iss_proofs_PartII.md` 之上继续。
> 目标：用 **AliSim + IQ-TREE 3.0.1** 重做 `manuscript.pdf` 中值得保留的实验，
> 并把每一项都接到解析闭式上。

---

## 0. 已定的决策（来自本轮确认）

| # | 决策 | 值 |
|---|---|---|
| D1 | 模拟引擎 | **AliSim 是唯一的模拟路径**（IQ-TREE 3.0.1 内置，含 HPC 模式）。不自己写模拟器 —— 可复现性要求用现成工具。AliSim 可控制 Q 矩阵、频率(`+F`)、自定义起始树(`-t`)。自建代码只负责：生成输入（树 + 设计密度给出的 s）与解析预测 |
| D2 | 推断引擎 | **IQ-TREE 3.0.1** 同一二进制 |
| D3 | 分类单元 | **N=4 为主体**；仅实验 2 额外加 N=8 / N=16 的**点迹图** |
| D4 | Item 1 的独立性 | **不预设结论**：出 (TL × L) 二维图，让数据说话 |
| D5 | Item 4 的功效 | **不做任何功效指标**，只诚实展示各统计量随单参数的变化趋势 |
| D6 | 采样 | **uniform-Iss 逆向采样**（Item 0 的设计密度），不靠 wall-clock 硬铺网格 |
| D7 | 运行方式 | 所有模拟**挂后台**；前台留给数学讨论 |
| D8 | 口径 | 主口径 = **Xia 口径 A**（$H_{FSS}=\mathbb E_{\mathrm{Mult}}[h]$，底数 $\log_2$）；所有图表标注口径，实验 2 同时给 B |
| D9 | Γ 离散化 | IQ-TREE 默认（等概率**均值**，Yang 1994）。另跑 `-gmedian` 对照臂，用于核对 Part I Table S3 |

---

## 1. 仓库结构

```
Iss_rebench/
├── PLAN.md                     本文件
├── README.md
├── src/issbench/
│   ├── entropy.py              h、H_FSS 口径 A/B、底数可切、精确多项和
│   ├── models.py               JC69/K2P/F81/F84/HKY85/TN93/GTR 的 Q 与谱
│   ├── spectral.py             Lemma 9/10 的 τ_k, w_k, ξ_k, Δ（Part II）
│   ├── gamma.py                Yang 离散 Γ（mean / median）+ 连续拉普拉斯
│   ├── trees.py                拓扑 + 枝长模式；zone transfer 路径；Newick 写出
│   ├── exact.py                Theorem 2 精确 E[I_ss]；N=4 的 256-pattern 精确 μ/σ²
│   ├── design.py               Item 0 的设计密度（闭式 + 精确数值逆）
│   ├── alisim.py               AliSim 调用封装（Newick → 比对，seed 管理，批处理）
│   ├── stats.py                Iss / p̄ / C-value / λ
│   ├── infer.py                IQ-TREE 调用 + P_true 判据（拓扑完全相同）
│   └── io.py                   结果落盘、断点续跑
├── scripts/
│   ├── s00_gate.py             ★ G0 门：AliSim 管线自检（必须通过才往下走）
│   ├── s01_design.py           Item 0
│   ├── s02_Lmap.py             Item 1
│   ├── s03_iss_pdist.py        Item 2
│   ├── s04_zone.py             Item 3
│   ├── s05_perturb.py          Item 4
│   └── s06_figures.py          Item 5
├── conf/*.yaml                 每个实验一套参数
├── results/                    单元结果 CSV + summary；大文件 gitignored
├── figs/                       图版
├── docs/                       Part I / Part II / 审计报告
│   └── audit/                  历史核验脚本与输出（冻结）
└── refs/                       refs_iss.bib / refs_iss.ris
```

---

## 2. Item 0 — uniform-Iss 设计密度

### 2.1 闭式推导

固定拓扑与枝长**模式** $\mathbf d^0$，整体缩放 $s$：$d_{ij}=s\,d^0_{ij}$。记

$$\nu(s)\;:=\;\frac{1}{\binom N2}\sum_{i<j}\mathbb E_r\!\Big[e^{-\frac43\,r\,s\,d^0_{ij}}\Big]$$

| 情形 | $\nu(s)$ | $\nu'(s)$ |
|---|---|---|
| 速率齐次 | $\frac1{\binom N2}\sum_{i<j}e^{-\frac43 s d^0_{ij}}$ | $-\frac{4}{3\binom N2}\sum_{i<j}d^0_{ij}e^{-\frac43 s d^0_{ij}}$ |
| $\Gamma(\alpha)$ | $\frac1{\binom N2}\sum_{i<j}\big(1+\frac{4 s d^0_{ij}}{3\alpha}\big)^{-\alpha}$ | $-\frac{4}{3\binom N2}\sum_{i<j}d^0_{ij}\big(1+\frac{4 s d^0_{ij}}{3\alpha}\big)^{-\alpha-1}$ |

**口径 A**（Xia）下 $\mu(s)=\mathbb E[I_{ss}]$ 的二阶闭式：

$$\mu(s)=1-\frac{3(N-1)}{2N}\cdot\frac{\nu(s)}{\,h(\boldsymbol\pi)-\frac{K-1}{2N}\,},\qquad
\mu'(s)=-\frac{3(N-1)}{2N\left(h(\boldsymbol\pi)-\frac{K-1}{2N}\right)}\,\nu'(s)$$

**设计密度**（Theorem 11）：

$$p_s(s)=\frac{\mu'(s)}{Z},\qquad Z=\mu(s_{\max})-\mu(s_{\min}),\qquad F(s)=\frac{\mu(s)-\mu_{\min}}{Z}$$

**逆采样**：抽 $U\sim\mathrm{Unif}[0,1]$，置 $\mu^*=\mu_{\min}+UZ$，解 $\mu(s)=\mu^*$。等价地先解 $\nu$：

$$\nu^*=\frac{2N\left(h(\boldsymbol\pi)-\frac{K-1}{2N}\right)(1-\mu^*)}{3(N-1)}$$

| 情形 | $s$ 的解析逆 |
|---|---|
| 齐次 | $s=\dfrac{3}{4\bar d^0}\ln\dfrac{1}{\nu^*}$，$\bar d^0=\frac1{\binom N2}\sum_{i<j}d^0_{ij}$ |
| Γ + 星树（所有 $d^0_{ij}$ 相同） | $s=\dfrac{3\alpha}{4d^0}\Big[(\nu^*)^{-1/\alpha}-1\Big]$ |
| Γ + 一般树 | Newton：$s\leftarrow s-(\mu(s)-\mu^*)/\mu'(s)$，解析导数，3–5 步 |

**无需 pilot**：归一化常数就是 $\mu$ 的增量本身。**不需要估计任何参数** ——
输入只有 $(\mathbf d^0,\ N,\ \alpha,\ \boldsymbol\pi)$，全部由构造成立。

**口径不变性**：$\mu_A=c\,\mu_B$（$c$ 为常数），故 $p_s^A\equiv p_s^B$ —— **设计密度本身与口径无关**；
只有目标 Iss 区间 $[\mu_{\min},\mu_{\max}]$ 需要按口径换算。

### 2.2 一个必须处理的问题：闭式 μ 是二阶渐近式

$|\mu_{\text{2阶}}-\mu_{\text{真}}|$ 在饱和区 $\approx5\times10^{-3}$，低分歧处可达 $0.06$
（Part I Table S2）。若直接用二阶 $\mu$ 做设计，Iss 边缘会**系统性倾斜**。

**方案（两用）**：
1. **operational**：N=4 用 `exact.py` 精确算 $\mu(s)$（256 pattern，秒级），数值反演 → 边缘严格平；
2. **诊断图**：把"用二阶 $\mu^{-1}$ 反解后 Iss 边缘的倾斜量"单独画成一张图 ——
   **这张图同时是采样器的误差条和 Theorem 1 有效区间的检验**（Part I Remark 3）。

### 2.3 产出

- `fig_design_density.png`：$\nu(s)$、$\mu(s)$、$p_s(s)$、$F(s)$ 四联图
- `fig_marginal_flatness.png`：精确 μ 反演 vs 二阶 μ 反演的 Iss 边缘对照
- `design.json`：每个实验实际使用的 $s$ 网格

**Commit**: `feat(design): closed-form uniform-Iss sampling density with exact-inversion fallback`

---

## 3. Item 1 — Iss 的位置与尺度：二维图

**不预设结论**，同时输出四张 (TL × L) 热图：**中位数 / 均值 / SD / CI 宽度**。

- 拓扑：4-taxon 对称 + 4-taxon 不对称（复刻 `manuscript.pdf` Fig 2A 的两列）
- 网格：$\mathrm{TL}\in[2,200]$（对数分布，约 24 档）× $L\in\{100,200,500,10^3,2\times10^3,5\times10^3,10^4\}$
- 每格 $n_{reps}$ 由**先验成本表**定：$n_{reps}=C(s)/(L\cdot \mathrm{SE}^2)$，$C(s)$ 从精确分布算出
- **N=8 / N=16 点迹图**（D3）：只取若干代表性的 (TL, L) 点叠加，验证关系是否跨 N 保持

**关键交叉检验**：位点独立 $\Rightarrow$ $n_{reps}\times L \equiv 1\times n_{reps}L$（Part II A5），
所以"跨重复的 SD"必须等于"同 TL 下不同 L 的 SD 按 $L^{-1/2}$ 外推"，两条路互证。

**Commit**: `feat(exp1): 2-D maps of Iss location and scale over (tree length x sites)`

---

## 4. Item 2 — JC69：Iss vs p̄ 与闭式对照

- 用 Item 0 的均匀-Iss 设计取 $s$（固定拓扑，扫 $s$），每点 $n_{reps}$ 次 AliSim
- 三条曲线叠：**MC 散点** / **Theorem 2 精确闭式**（实线）/ **Theorem 1 二阶仿射**（虚线，展示失效区）
- 口径：主线用 **A**；副图用 **B**，两者比值 $c$ 实测标出
- **N=8 / N=16 点迹图**（D3，用户指定）：同图叠上，展示关系跨 taxa 保持
- 附：$\bar p$ 的解析 $\frac34(1-\nu)$ 叠加

**Commit**: `feat(exp2): Iss vs p-distance with exact and second-order closed forms (N=4/8/16)`

---

## 5. Item 3 — Felsenstein → Farris 转移

**路径**（用户指定，已证明与 Part II Theorem 9 的路径**逐点等价**）：
$$(a,b,c,d)=(L-t,\;S,\;L,\;S+t),\qquad t:0\to L-S$$
- $t=0$：$(L,S,L,S)$ —— 长枝 A、C 被内枝 $r$ 隔开 ⇒ **Felsenstein zone**
- $t=L-S$：$(S,S,L,L)$ —— 长枝 C、D 相邻 ⇒ **Farris zone**

**与 Part II 的对接（已核）**：两参数化在任一 $t$ 给出的 pairwise 距离多重集完全相同
（差一个 clade 内的交换），而 $H$ 对叶标签对称、两个 clade 可交换
$\Rightarrow$ **精确 $\mathbb E[H]$ 沿两条路径逐点相同**。故 Theorem 9 的三条结论逐字适用：

1. $\sum_{i<j}d_{ij}=3(a+b+c+d)+4r=6L+6S+4r$ 恒定
2. $\bar p$ 与 $\mathbb E[I_{ss}]$ **不**恒定（Remark 13）
3. $r=0$ 时两端多重集 $\{2L,2S,(L+S)\times4\}$ 完全重合 ⇒ $I_{ss}$ 严格相等
4. $O(r)$ 间隙 $=\frac{N-1}{4Nh}\cdot\frac{4r}{3}(e^{-4L/3}-e^{-4S/3})^2$

**每次运行**：沿 $t$ 取 ~17 点 × 每点 $n_{reps}$ 次 AliSim，得
(i) Iss（实测 + 精确闭式）；(ii) $P_{true}$，**同时跑 ML（IQ-TREE）与 MP、NJ**（Part II Remark 14 要求两种都报）；
(iii) 四点法判据（用于与 Part I §S9 的表对话 —— 见 §9 风险 R2）。

**主图**（$x=$Iss, $y=P_{true}$，点按 $t$ 着色）：同一 Iss 值上 $P_{true}$ 从 ≈1 散到 ≈0。
**副图**：$\mathbb E[I_{ss}]$ vs $t$ 与实测值叠（检验 Theorem 9.2 与 9.4）。

**Commit**: `feat(exp3): Felsenstein->Farris transfer; Iss flat while P_true spans the range`

---

## 6. Item 4 — 单参数扰动 vs C-value / λ

**不做功效指标**（D5）。一次只动一个参数，画各统计量的趋势线与置信带。

| 旋钮 | 档位 | 固定其余 |
|---|---|---|
| $\alpha$（Γ 形状） | $\infty,4,2,1,0.5,0.25$ | JC69、星树/对称树、TL 固定 |
| $\kappa$（ts/tv） | $1,2,5,10,20$ | HKY、均匀 $\pi$、TL 固定 |
| $\pi$ 偏斜 | 均匀 → $(.35,.15,.15,.35)$ → $(.7,.1,.1,.1)$ | HKY $\kappa{=}5$、TL 固定 |
| 拓扑不平衡度 | 对称 → 中间 → pectinate | JC69、TL 固定 |
| 内枝 $r$ | 0.005 → 0.5 | JC69、pendant 固定 |
| $L$ | $10^2\to10^4$ | 其余固定 |

每个统计量：**Iss（口径 A）**、**$\bar p$**、**C-value**、**$\lambda$**。

⚠ **先决条件**：C-value 与 λ 的定义必须先钉死并写进 `stats.py` 的 docstring（附出处 DOI）：
- C-value：Struck et al. 2008, MPE 48:628–645 `10.1016/j.ympev.2008.05.015`
- λ：Jermiin & Misof 2020, bioRxiv `10.1101/2020.01.03.897097`（**仅预印本**）
- 两者都写进 `refs/refs_iss.bib`（已在 B 组）

**Commit**: `feat(exp4): one-at-a-time perturbation trends for Iss / C-value / lambda`

---

## 7. Item 5 — 图版

按 `ref2.png` 的风格：四面板、衬线小标题（"A. …"）、无衬线轴标签、低饱和配色
（黑 / 深红 / 橙 / 蓝 / 绿）、细边框、浅灰网格、虚线参考线、米色标定区。

**Commit**: `feat(figs): ref-style plate for the simulation results`

---

## 8. 运行与调度

- 所有模拟经 `run_in_background` 启动，日志落 `logs/<exp>/<cell>.log`
- 每个 (cell) 独立可重跑；结果存在即跳过（断点续跑）
- **seed 纪律**：`--seed` = 确定性函数 `f(exp, cell, rep)`，写进结果表；任何重跑 bit-identical
- 预算控制：$P_{true}$ 只在过渡区（$P\in[0.1,0.99]$）密铺；区间外答案已知

**AliSim 调用骨架**（待 G0 验证）：
```bash
iqtree3 --alisim <prefix> -t <tree.nwk> -m "JC69+G4{1.0}" \
        --length <L> --num-alignments <n> --seed <s> -af fasta
```

---

## 9. 阶段门与风险

### G0 门（**必须先过**）
`s00_gate.py`：建 JC69 星树 N=16 的 Newick，跑 AliSim，检验
1. `--seed` 可复现（两次运行 bit-identical）
2. `--branch-scale` 确为**乘法**（否则改用直接写 Newick）
3. 实测 $\bar p$ 与解析 $\frac34(1-e^{-8d/3})$ 一致
4. 实测 Iss 与 **Part I Theorem 2 精确闭式**一致（口径 A 与 B 各算一遍）
5. Γ 离散化确认为**均值**方案（与 `-gmedian` 的结果可区分）

### 风险

| # | 风险 | 处置 |
|---|---|---|
| R1 | AliSim 的 `+G4` 未按预期离散化 | G0 第 5 项钉死；必要时改用 `--site-rate MODEL` 或直接喂连续 Γ |
| R2 | **Part I §S9 的 $P_{true}$ 表标注为 "four-point condition"，但独立复算显示四点法在 $TL{=}32,\alpha{=}0.5$ 给 0.895，而表里是 1.000（ML 给 1.000）** | 实验 3 用 AliSim + IQ-TREE 重跑该表，同时保留四点法。**这是数学讨论的开放项** |
| R3 | uniform-Iss 设计在低分歧处依赖二阶 μ | §2.2 的双轨方案 |
| R4 | $P_{true}$ 的 wall-clock 在过渡区放大 | 先验成本表 + 只在过渡区密铺 |
| R5 | C-value / λ 定义歧义 | 实现前先对照原文，写死 docstring |

---

## 9b. 实验 6（新增，来自 gpt6 审查的拾遗）

**条件独立检验 —— 项目原先的 $R^2$ 论证回答不了这个问题。**

在控制平均 $p$-距离之后，$I_{ss}$ 是否仍携带独立预测信息？

$$\Pr(Y\mid ar p, I_{ss})\quad	ext{vs}\quad\Pr(Y\mid ar p)$$

$Y$ 取两个真值：(i) 简约信息位点比例；(ii) 拓扑恢复成功与否。

设计：复用实验 2 的 uniform-Iss 采样点（同一 $ar p$ 上有多个 $I_{ss}$ 的条件由设计密度保证），
在每个 $ar p$ 分箱内比较两个 logistic 模型的拟合与交叉验证 AUC。

**为什么必须做**：项目的"冗余"主张建立在二阶论证上（$\E[I_{ss}]$ 只依赖 $\langle P_{jj}angle$），
而**二阶之上** $\E[I_{ss}]$ 依赖位点模式分布的三阶、四阶矩，那些不由成对同一性决定。
$N{=}4$ 时 2:2 是简约信息位点、3:1 不是——两种模式熵不同。所以 $I_{ss}$ 原则上能读到
$ar p$ 读不到的模式结构。控制 $ar p$ 之后条件独立与否，是"冗余"的唯一判定点。

---

## 10. 提交序列总览

```
ab7df45  chore: snapshot pre-existing audit material and initialize repository
8efe887  chore(repo): consolidate audit scripts into docs/audit, ...
         ← 以下为本次计划
         docs: add PLAN.md
         feat(io)+feat(models)+feat(spectral): package skeleton
         test: AliSim gate G0                       ← 门
         feat(design): uniform-Iss sampling density ← Item 0
         feat(exp1): 2-D maps of Iss location and scale
         feat(exp2): Iss vs p-distance with closed forms
         feat(exp3): Felsenstein->Farris transfer
         feat(exp4): perturbation trends vs C-value / lambda
         feat(exp5): conditional-independence test (Pr(Y|pbar,Iss) vs Pr(Y|pbar))
         feat(figs): ref-style plate
```
