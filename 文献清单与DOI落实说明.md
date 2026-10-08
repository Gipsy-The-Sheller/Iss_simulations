# 文献清单与 DOI 落实说明

> 覆盖 `Iss方法学批判_技术要点汇总..md` 与 `Iss饱和指数_方法学批判与复现记录.md` 中**点名的全部文献**。
> 输出：`refs_iss.bib`（BibTeX，39 条）与 `refs_iss.ris`（RIS，39 条）。
> 全部元数据由 `verify/make_refs.py` 从 **Crossref API 实时抓取**（缓存于 `verify/crossref_cache.json`），
> 无一条凭记忆录入。每条均带 `note` 字段标注它在哪一节被提到、起什么作用。
> 抓取日期：2026-10-08。

---

## 1. 清单总览（39 条）

| 组 | 条数 | 内容 |
|---|---|---|
| **A 核心对象：Xia / DAMBE** | 8 | 批判对象原文、DAMBE 四代、starless bias、教科书、Kubatko 书评、PhyPA |
| **B 替代与对照诊断量** | 9 | C 值、nRCFV、λ、SatuTe、PhyloMAd（论文+软件）、TreSpEx、BaCoCa、PI |
| **C 先例与相关综述** | 6 | Verbruggen 2007、Letsch & Kjer 2011、a priori 回顾、Robust Phylogenetics、MBE 40 年、clawRxiv |
| **D 软件与工具** | 10 | IQ-TREE 2、AliSim、PAML 4、PhyloBayes 3、BEAST 1.10、TreeTime、Nextstrain、PhyML 3.0、MAFFT 7、MUSCLE |
| **E 病毒定年（TRAD）** | 4 | Xia 2021/2022/2023/2025 四代 TRAD |
| **F 补充（核验新增）** | 2 | Xia & Lemey 2009（Iss.c 规范引用）、Yang 1994（离散 Γ 原始文献） |

---

## 2. ⚠ 需要修订 md 的地方（落实过程中发现的 5 处）

### 2.1 「Verbruggen et al. 2007」的作者序与实际文献

汇总 §3.1 / §11.7 与记录 §2.1 均写「**Verbruggen et al. 2007（Mol Phylogenet Evol）**」。

实际对应的是：
> **Verbruggen H**, Leliaert F, Maggs CA, Shimada S, Schils T, Provan J, Booth D, Murphy S, De Clerck O, Littler DS, Littler MM, Coppejans E (2007)
> *Species boundaries and phylogenetic relationships within the green algal genus Codium (Bryopsidales) based on plastid DNA sequences.*
> Mol Phylogenet Evol **44**(1):240–254. `10.1016/j.ympev.2007.01.009`

- ✅ **DOI 落实、作者序正确**（Verbruggen 确为第一作者）。
- ⚠ **但"该文指出 Iss.c 依赖 taxa 数与序列长度、因而是数据集专属无法制表"这一具体论断，我未能核实。**
  检索只确认该文使用了 DAMBE 的饱和检验并报告 Iss < Iss.c。**必须调出全文核对原句**，否则在正文里
  以"Verbruggen et al. 2007 已指出"引用是引用失实。
- 注意区分：另有 **Leliaert F, De Clerck O, Verbruggen H, … (2007)** *Molecular phylogeny of the
  Siphonocladales…* MPE 44(3):1237–1256 `10.1016/j.ympev.2007.04.016`，Verbruggen 是**第三作者**，
  常被误引为 "Leliaert et al. 2007"。**若引文要写 "et al."，只有前一篇才是 Verbruggen et al.**

### 2.2 「SatuRation / nRCFV_Reader」（汇总 §10.2）

- ✅ **nRCFV_Reader 真实存在**（RCFV 与 nRCFV 的计算工具，github.com/JFFleming/RCFV_Reader），
  对应文献 = Fleming & Struck 2023（已收入 B 组）。
- ❌ **"SatuRation" 未检索到任何对应工具或论文。** 检索到的近似项是 **SatuTe**（同组已收），
  二者字母相近，疑为笔误或记忆混淆。**建议改名为 SatuTe，或删去该词。**

### 2.3 「Jermiin & Misof 2020」（汇总 §8.2、§11.7）

- ✅ 文献存在：*Measuring Historical and Compositional Signals in Phylogenetic Data*，`10.1101/2020.01.03.894097`
- ⚠ **只有 bioRxiv 预印本，未检索到期刊版本。** 汇总 §11.7 把它与正式论文并列在"文献"表里，
  正文引用时应标注 "(preprint)"，否则会被认为引用了不存在的东西。λ 的定义（可变位点比例 ÷ 无限时间
  期望比例）与 δCFS 均出自此文。

### 2.4 「clawRxiv 2604.01176 (2026)」（汇总 §11.7、记录 §9）

- ✅ **确实存在**：*The Substitution Saturation Threshold: Phylogenetic Signal Becomes Unrecoverable
  Beyond 0.8 Substitutions Per Site for Protein-Coding Genes*，https://www.clawrxiv.io/abs/2604.01176
- ❌ **无 DOI**，不在 Crossref 中。
- 🔎 **新证据**：clawRxiv 同一目录下还有 *"Research Agent: An Executable Skill for Autonomous Web
  Research"*、*"TAN-POLARITY v4: A Pre-Validation Framework…"* 等条目 —— 该平台是 **AI-agent 生成的
  可执行预印本服务器**，不是学术期刊体系的一部分。
  **这为记录 §9"未采信其数字"的判断提供了直接支撑**，可直接写进正文脚注。

### 2.5 「Algorithms for Molecular Biology（a priori 数据质量评估综述）」（汇总 §9.4）

- ✅ 落实为：**Misof B, Meusemann K, von Reumont BM, Kück P, Prohaska SJ, Stadler PF (2014)**
  *A priori assessment of data quality in molecular phylogenetics.* Algorithms for Molecular Biology **9**:22.
  `10.1186/s13015-014-0022-4`
- 年份是 **2014** 不是 2021，汇总未写年份（只写了刊名），补上即可。

### 2.6 其他已确认的年份/卷页更正

| md 写法 | Crossref 实际 | 说明 |
|---|---|---|
| Duchêne et al. **2021** (Syst Biol) | Syst Biol **71**(3):676–689；online 2021、issue 2022 | 引 2021 可，但建议写 `2021` 并给 71(3) |
| Kubatko 书评 | Syst Biol **70**(2):408–411, **2020** | Crossref 标 2020（issue 2021） |
| Duchêne et al. 2021 提到的 PhyloMAd | **Duchêne, Duchêne & Ho 2018**, Bioinformatics 34(13):2300–2301 | 软件论文与 2021 论文是两篇，别混 |
| Xia 2020《A Mathematical Primer…》 | 出版社 **Apple Academic Press**，`10.1201/9780429425875` | CRC/T&F 代理发行 |

---

## 3. 我额外补进去的两条（F 组，核验报告承重）

| key | 文献 | 为什么必须收 |
|---|---|---|
| `xia2009handbook` | Xia X & Lemey P (2009) *Assessing substitution saturation with DAMBE*. In: The Phylogenetic Handbook, CUP, pp. 615–630. `10.1017/CBO9780511819049.022` | **Iss.c 被广泛使用的规范引用**（比 2003 原文更多被实践者引）。两份 md 通篇未引，但审稿人一定会问"你用哪个版本的检验"。 |
| `yang1994discretegamma` | Yang Z (1994) *Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods.* J Mol Evol **39**(3):306–314. `10.1007/BF00160154` | **离散 Γ 的标准等概率均值离散化出处**。核验报告 A1 指出笔记用的是中位数方案，与本文冲突时**必须引此文**界定"标准 Γ4"是什么。 |

---

## 4. 未纳入的两类

1. **Xia 2003 原文自己的 40 条参考文献**（Archie 1989、Faith 1991、Swofford 1993、Steel et al. 1993/1995、
   Kuhner & Felsenstein 1994、Yang 1996a/b/2000、Ritland & Clegg 1990 等）—— 两份 md 均未提及，
   按"md 中点名"的范围未收。若论文要复述 Xia 的"四类方法"综述，需要另行补齐，可再跑一次脚本。
2. **纯软件名**（SEAVIEW、TempEst、nRCFV 的多个下游工具）—— 只在行文中一笔带过，未点名文献。

---

## 5. 复现方式

```bash
python verify/make_refs.py        # 从 Crossref 重新抓取并生成 refs_iss.bib / refs_iss.ris
```

- 元数据缓存在 `verify/crossref_cache.json`；删掉该文件即强制重新联网抓取。
- 需要改动（加条目、改正条目、换 key）时**改 `make_refs.py` 里 ENTRIES 与 NOTE 两张表**，不要手改输出文件。
- 每个条目的 `note` 字段已写明它在 md 的哪一节出现，方便正文写作时回溯。
- BibTeX 的 `note` 是中文，RIS 的 `N1` 同；**投稿前记得清空这两个字段**，否则会跟着参考文献一起印出来。
