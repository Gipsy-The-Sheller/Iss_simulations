# Research Square 预印本包

交付物是**一个 PDF**。Research Square 只接受 `.doc` / `.docx` / `.pdf`（编辑政策 §3），LaTeX 源文件不作为提交件，所以这个目录不需要像 `MPE_sumb/` 那样自包含。

## 目录内容

| 文件 | 说明 |
|---|---|
| `manuscript.tex` | 预印本正文（`article` 11pt，保留行号），由 `../manuscript.tex` 生成 |
| `manuscript.pdf` | 上传用的成品，35 页（原稿 34 页，多出的是新增声明） |
| `manuscript.bbl` | BibTeX 输出，按仓库惯例保留 |
| `_abstract.txt` | 纯文本摘要（328 词，三段），提交表单里粘贴用 |
| `_build_rs.py` `_verify.py` | 生成与核对脚本，投稿时不必上传 |
| `../figures/` `../refs/` | **不复制**：正文通过 `\graphicspath{{../}}` 与 `../refs/` 指回仓库根目录 |

## 相对原稿（`../manuscript.tex`）改了什么

只有六处，正文一字未动——`_verify.py` 会逐字核对这一点。

1. **删起草注释**：文件头 25 行（结构备注与标题候选）、9 条 `% ====` 分隔线、文末 14 行 `NOTE (do not reintroduce)` 块。原稿里**没有行尾注释**，行内出现的 37 个 `%` 全是转义的 `\%`（百分号），不受影响。
2. **两处路径**：导言区加 `\graphicspath{{../}}`，正文里 `figures/xxx.png` 的写法不动；`\bibliography{refs/...}` 改为 `../refs/...`。
3. **作者脚注**补 `Correspondence: zhi-jie.xu@dhu.edu.cn`。
4. **摘要第二句加一句新工作定位语**："We characterize the index exactly, and test that characterization in new simulations and in re-analyses of two published data sets."（正文的 `_abstract.txt`、表单里的摘要，与 PDF 一致。）理由见下面「为什么加这句」。**标题未动**——它是这篇稿子的引用身份，要与 MPE 投稿版一致。
5. **`Use of AI-assisted tools`** 小节，加在 Material and methods 末尾（Research Square 编辑政策 §1.3 要求 LLM 使用须在 Methods 中说明）。措辞为：Claude Code (Anthropic) 用于文稿的撰写与编辑，以及数学结果的推导与检验；作者对内容负责，工具不作为作者。
6. **三段声明**（Data availability 之后）：CRediT、Declaration of competing interest、Funding。措辞逐字取自 `MPE_sumb/_build_mpe.py`，两版不会各说各话。

**没有改的**：章节顺序仍沿用原稿（Results 在 Material and methods 之前——政策 §1.2 明确写 "alternative formats are also accepted"）；保留行号；拼写仍是英式；参考文献样式仍是 `plainnat`。

## 提交前需要你确认

- [ ] **AI 声明措辞**（§1.3 强制，但要你签字）：现按"写作 + 数学证明的推导与检验"写，未提代码实现。要加点名到具体定理、或要收窄/放宽范围，改 `_build_rs.py` 里的 `AI_DISCLOSURE` 后重跑即可。
- [ ] **利益冲突**：PDF 里已声明"无"，**表单里还有一处勾选/填写**，两处要一致。
- [ ] **Funding**：现写"无专项资助"，有基金请替换。
- [ ] **通讯作者邮箱**：占位 `zhi-jie.xu@dhu.edu.cn`（与 MPE 版同一处待确认）。
- [ ] **时点**：政策 §6 明确"已在期刊发表的稿件不会被接收"。MPE 正在投，挂预印本会不会影响那边，由你判断（Elsevier 一般允许预印本）。
- [ ] **表单内容与 PDF 一致**：题名、摘要用 `_abstract.txt`，关键词用原稿 `Key words` 那 8 个——表单里填的内容**会原样公开**。摘要 328 词，RS 没有明文上限；若表单卡长度，`MPE_sumb/_abstract.tex`（248 词）是现成的备选。
- [ ] 若要单独上传 figure files（表单可选），用 `../figures/` 里的 12 张 PNG，均 ≥300 dpi。
- [ ] Article type 选 **Research Article**；学科领域按你的判断选（§1.2 要求完整的 Methods 与 Results，本文具备）。

## 重新编译

```bash
cd RS_preprint
pdflatex manuscript && bibtex manuscript && pdflatex manuscript && pdflatex manuscript
# 或者：python _build_rs.py   （从 ../manuscript.tex 重新生成 .tex 与 _abstract.txt）
python _verify.py             （核对生成物是否最新、正文是否只有上述改动）
```

## 为什么加那句定位语

RS 编辑政策 §1.2 只收 "research articles with complete methods and results sections"，明确排除 "Theories, commentaries, and non-systematic reviews"——而这篇稿子的形状（对一篇已发表论文的方法学批评）正是筛查编辑可能误归成 commentary 的那种。**bioRxiv 那条更硬的规则就是前车之鉴**：

> bioRxiv is intended for the rapid sharing of new research and not for commentary that simply contests or rebuts research published in journals. However, articles that augment or challenge earlier work with bona fide new data and/or complete re-analyses of published data sets may be deemed appropriate for the Contradictory Results section.

稿子其实完全符合它的豁免条款（新定理 + 108,000 条新模拟 + 对已发表数据集的完整重分析），只是没把自己呈现成那样。摘要第二句现在明说了这件事，标题保持不动。若 RS 也卡在体裁上，退路是 EcoEvoRxiv——那边明确欢迎 opinion/comment 类（政策：opinion, perspective, and comment pieces），不动一个字就能挂；但注意它**不接受已在别的平台发布的稿子**，两边只能选一个。

## 两点提醒

- 政策 §6 的筛查会拒掉缺图、缺参考文献列表的稿子，也会对"结论很强但数据不可得"的稿子多看一眼——文内的 Data availability（Zenodo DOI + GitHub）正好是这一条的解药，别把它删了。
- 预印本一经发布即获 DOI，不可撤除（§8/§13）；CC-BY 4.0 自动适用（§14）。
