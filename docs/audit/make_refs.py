# -*- coding: utf-8 -*-
"""Build refs_iss.bib + refs_iss.ris for every work named in the two .md files.
Metadata is pulled live from Crossref (api.crossref.org/works/{doi}) so nothing
is typed from memory."""
import os, json, urllib.request, urllib.parse, time, re, html

UA = {"User-Agent": "IssCritiqueBib/1.0 (mailto:researcher@example.org)"}
NOTE = {  # where the .md files mention it / why it is in the list
 "xia2003iss": "批判对象原文。汇总 §1.2、记录 §1",
 "xia2001dambe": "DAMBE 初版；记录 §8 被引量",
 "xia2013dambe5": "DAMBE5；记录 §8",
 "xia2018dambe7": "DAMBE7；记录 §8",
 "xia2018starless": "starless bias；汇总 §3.6、记录 §3.3（α 与拓扑强混淆）",
 "xia2020primer": "Xia 教科书；记录 §8",
 "kubatko2021review": "Kubatko 书评；记录 §8",
 "xia2016phypa": "PhyPA；记录 §8（无第三方评测）",
 "struck2008cfactor": "C 值原始文献；汇总 §8.1",
 "kuck2014bacoca": "BaCoCa 实现；汇总 §8.1、§10.2",
 "struck2014trespex": "TreSpEx；汇总 §10.2",
 "jermiin2020signal": "λ 历史信号/组成距离；汇总 §8.2、§11.7（注：仅 bioRxiv 预印本）",
 "fleming2023nrcfv": "nRCFV；汇总 §11.7、记录 §2.3",
 "manuel2025satute": "SatuTe；汇总 §9.1（分支级饱和检验，本文正面落点）",
 "duchene2021saturation": "汇总 §10.3、§11.7 —— 对 √N 与 L 回归，模拟按'长树长'操作化饱和",
 "duchene2018phylomad": "PhyloMAd 软件论文；汇总 §7、§9.1",
 "townsend2007pi": "phylogenetic informativeness；记录 §7 替代方案表",
 "verbruggen2007codium": "汇总 §3.1、§11.7 —— Iss.c 依赖 taxa 数与 L、无法制表",
 "letsch2011rrna": "汇总 §11.7 —— 饱和非位点固有属性，取决于系统发育层级",
 "amb2014apriori": "汇总 §9.4 —— 同塑性检测不得预设待推导的树",
 "liu2025robustphylo": "汇总 §9.4 —— 树依赖方法的循环性",
 "clawrxiv260401176": "汇总 §11.7、记录 §9 —— 非主流预印本平台，无 DOI",
 "minh2020iqtree2": "汇总 §11.8 推断引擎",
 "lytrong2022alisim": "汇总 §11.8 模拟引擎（确定性种子）",
 "yang2007paml4": "汇总 §1.4 EVOLVER / PAML",
 "lartillot2009phylobayes": "记录 §7 后验预测检验",
 "suchard2018beast": "记录 §8 病毒定年生态",
 "sagulenko2018treetime": "记录 §8 病毒定年生态",
 "hadfield2018nextstrain": "记录 §8 病毒定年生态",
 "guindon2010phyml": "汇总 §8 记录 §8 —— DAMBE 建树模块封装对象",
 "katoh2013mafft": "同上（比对）",
 "edgar2004muscle": "同上（比对）",
 "xia2021trad": "TRAD 原始（恒定速率）；记录 §8",
 "xia2022tradimproved": "TRAD 改进（速率随时间线性）；记录 §8",
 "xia2023trad": "TRAD 应用；记录 §8",
 "xia2025tradhostswitch": "TRAD 最新（S 形速率，宿主转换）",
 "russo2024forty": "记录 §8 —— DAMBE 与 PAML/IQ-TREE/SEAVIEW 并列的出处",
 "xia2009handbook": "【补充】Iss.c 的规范引用（The Phylogenetic Handbook 章），两份 md 均未引",
 "yang1994discretegamma": "【补充】离散 Γ 的标准等概率均值离散化出处 —— 见核验报告 A1",
}

ENTRIES = [
 # ---- A 核心对象 ----
 ("A 核心对象：Xia / DAMBE", "xia2003iss", "10.1016/S1055-7903(02)00326-3"),
 ("A 核心对象：Xia / DAMBE", "xia2001dambe", "10.1093/jhered/92.4.371"),
 ("A 核心对象：Xia / DAMBE", "xia2013dambe5", "10.1093/molbev/mst064"),
 ("A 核心对象：Xia / DAMBE", "xia2018dambe7", "10.1093/molbev/msy073"),
 ("A 核心对象：Xia / DAMBE", "xia2018starless", "10.3934/genet.2018.4.212"),
 ("A 核心对象：Xia / DAMBE", "xia2020primer", "10.1201/9780429425875"),
 ("A 核心对象：Xia / DAMBE", "kubatko2021review", "10.1093/sysbio/syaa082"),
 ("A 核心对象：Xia / DAMBE", "xia2016phypa", "10.1016/j.ympev.2016.07.001"),
 # ---- B 替代与对照诊断量 ----
 ("B 替代与对照诊断量", "struck2008cfactor", "10.1016/j.ympev.2008.05.015"),
 ("B 替代与对照诊断量", "kuck2014bacoca", "10.1016/j.ympev.2013.09.011"),
 ("B 替代与对照诊断量", "struck2014trespex", "10.4137/EBO.S14239"),
 ("B 替代与对照诊断量", "jermiin2020signal", "10.1101/2020.01.03.894097"),
 ("B 替代与对照诊断量", "fleming2023nrcfv", "10.1186/s12859-023-05270-8"),
 ("B 替代与对照诊断量", "manuel2025satute", "10.1093/molbev/msaf090"),
 ("B 替代与对照诊断量", "duchene2021saturation", "10.1093/sysbio/syab075"),
 ("B 替代与对照诊断量", "duchene2018phylomad", "10.1093/bioinformatics/bty103"),
 ("B 替代与对照诊断量", "townsend2007pi", "10.1080/10635150701311362"),
 # ---- C 先例与相关综述 ----
 ("C 先例与相关综述", "verbruggen2007codium", "10.1016/j.ympev.2007.01.009"),
 ("C 先例与相关综述", "letsch2011rrna", "10.1186/1471-2148-11-146"),
 ("C 先例与相关综述", "amb2014apriori", "10.1186/s13015-014-0022-4"),
 ("C 先例与相关综述", "liu2025robustphylo", "10.1101/2025.04.01.646540"),
 ("C 先例与相关综述", "clawrxiv260401176", None),
 ("C 先例与相关综述", "russo2024forty", "10.1093/molbev/msad264"),
 # ---- D 软件与工具 ----
 ("D 软件与工具", "minh2020iqtree2", "10.1093/molbev/msaa015"),
 ("D 软件与工具", "lytrong2022alisim", "10.1093/molbev/msac092"),
 ("D 软件与工具", "yang2007paml4", "10.1093/molbev/msm088"),
 ("D 软件与工具", "lartillot2009phylobayes", "10.1093/bioinformatics/btp368"),
 ("D 软件与工具", "suchard2018beast", "10.1093/ve/vey016"),
 ("D 软件与工具", "sagulenko2018treetime", "10.1093/ve/vex042"),
 ("D 软件与工具", "hadfield2018nextstrain", "10.1093/bioinformatics/bty407"),
 ("D 软件与工具", "guindon2010phyml", "10.1093/sysbio/syq010"),
 ("D 软件与工具", "katoh2013mafft", "10.1093/molbev/mst010"),
 ("D 软件与工具", "edgar2004muscle", "10.1093/nar/gkh340"),
 # ---- E 病毒定年（Xia 系）----
 ("E 病毒定年 (TRAD)", "xia2021trad", "10.3390/v13091790"),
 ("E 病毒定年 (TRAD)", "xia2022tradimproved", "10.1007/978-3-662-65902-1_19"),
 ("E 病毒定年 (TRAD)", "xia2023trad", "10.3390/v15030684"),
 ("E 病毒定年 (TRAD)", "xia2025tradhostswitch", "10.1093/gbe/evaf134"),
 # ---- F 补充 ----
 ("F 补充（核验新增，直接承重）", "xia2009handbook", "10.1017/CBO9780511819049.022"),
 ("F 补充（核验新增，直接承重）", "yang1994discretegamma", "10.1007/BF00160154"),
]

CACHE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "crossref_cache.json")
try:
    CACHE = json.load(open(CACHE_PATH, encoding="utf-8"))
except Exception:
    CACHE = {}

def meta(doi):
    if doi in CACHE:
        d = CACHE[doi]
    else:
        u = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
        last = None
        for attempt in range(6):
            try:
                d = json.loads(urllib.request.urlopen(
                    urllib.request.Request(u, headers=UA), timeout=60).read())["message"]
                break
            except Exception as e:
                last = e; time.sleep(3 * (attempt + 1))
        else:
            raise SystemExit("FAILED to fetch %s : %s" % (doi, last))
        CACHE[doi] = d
        json.dump(CACHE, open(CACHE_PATH, "w", encoding="utf-8"), ensure_ascii=False)
    def dec(s): return html.unescape(re.sub(r"<[^>]+>", "", s or ""))
    return dict(
        doi=doi, title=dec((d.get("title") or [""])[0]),
        subtype=d.get("type"), journal=dec((d.get("container-title") or [""])[0]) if d.get("container-title") else "",
        publisher=d.get("publisher", ""), volume=d.get("volume", ""), issue=d.get("issue", ""),
        page=d.get("page", "") or d.get("article-number", ""), year=(d.get("issued", {}).get("date-parts") or [[None]])[0][0],
        authors=[(a.get("family", ""), a.get("given", "")) for a in (d.get("author") or [])],
        editors=[(a.get("family", ""), a.get("given", "")) for a in (d.get("editor") or [])],
        isbn=d.get("ISBN") or [], url=d.get("URL", ""))

def bibtex(key, m, note):
    typ = {"journal-article": "article", "book": "book", "reference-book": "book",
           "edited-book": "book", "monograph": "book", "book-chapter": "incollection",
           "posted-content": "misc", "proceedings-article": "inproceedings"}.get(m["subtype"], "misc")
    L = ["@%s{%s," % (typ, key)]
    if m["authors"]:
        L.append("  author  = {%s}," % " and ".join("%s, %s" % (f, g) if g else f for f, g in m["authors"]))
    L.append("  title   = {%s}," % m["title"].rstrip("."))
    if typ == "article":
        L.append("  journal = {%s}," % m["journal"])
    elif typ == "incollection":
        L.append("  booktitle = {%s}," % m["journal"])
    if m["volume"]: L.append("  volume  = {%s}," % m["volume"])
    if m["issue"]:  L.append("  number  = {%s}," % m["issue"])
    if m["page"]:   L.append("  pages   = {%s}," % m["page"].replace("-", "--"))
    L.append("  year    = {%s}," % m["year"])
    if m["publisher"]: L.append("  publisher = {%s}," % m["publisher"])
    if m["isbn"]:   L.append("  isbn    = {%s}," % m["isbn"][0])
    L.append("  doi     = {%s}," % m["doi"])
    L.append("  note    = {%s}," % note)
    L.append("}")
    return "\n".join(L)

RIS_TY = {"journal-article": "JOUR", "book": "BOOK", "book-chapter": "CHAP",
          "posted-content": "UNPB", "proceedings-article": "CONF"}
def ris(key, m, note):
    L = ["TY  - %s" % RIS_TY.get(m["subtype"], "GEN")]
    for f, g in m["authors"]:
        L.append("AU  - %s" % (("%s, %s" % (f, g)) if g else f))
    L.append("TI  - %s" % m["title"])
    if m["subtype"] == "article":        L.append("T2  - %s" % m["journal"])
    elif m["subtype"] == "book-chapter": L.append("T2  - %s" % m["journal"])
    if m["volume"]: L.append("VL  - %s" % m["volume"])
    if m["issue"]:  L.append("IS  - %s" % m["issue"])
    if m["page"]:
        sp = m["page"].split("-")[0].strip()
        L.append("SP  - %s" % sp)
        if "-" in m["page"]: L.append("EP  - %s" % m["page"].split("-")[-1].strip())
    L.append("PY  - %s" % m["year"])
    L.append("PB  - %s" % m["publisher"])
    L.append("DO  - %s" % m["doi"])
    L.append("UR  - %s" % m["url"])
    L.append("N1  - %s  [%s]" % (note, key))
    L.append("ER  - ")
    return "\n".join(L)

bib, risout, cur = [], [], None
bib.append("% refs_iss.bib -- 由 verify/make_refs.py 从 Crossref 自动生成，请勿手改。\n"
           "% 覆盖两份 .md 中点名的全部文献。抓取日期 2026-10-08。\n")
risout.append("; refs_iss.ris -- 由 verify/make_refs.py 从 Crossref 自动生成，请勿手改。\n"
              "; 覆盖两份 .md 中点名的全部文献。抓取日期 2026-10-08。\n")
HDR = "%s ======================= %s =======================\n"

for group, key, doi in ENTRIES:
    if group != cur:
        cur = group
        bib.append("\n" + HDR % ("%%", group))
        risout.append("\n" + HDR % (";", group))
    note = NOTE.get(key, "")

    if doi is None:                       # the one DOI-less entry (clawRxiv)
        title = ("The Substitution Saturation Threshold: Phylogenetic Signal Becomes "
                 "Unrecoverable Beyond 0.8 Substitutions Per Site for Protein-Coding Genes")
        url = "https://www.clawrxiv.io/abs/2604.01176"
        bib.append("@misc{%s,\n"
                   "  author       = {{clawRxiv}},\n"
                   "  title        = {%s},\n"
                   "  year         = {2026},\n"
                   "  howpublished = {clawRxiv preprint 2604.01176},\n"
                   "  url          = {%s},\n"
                   "  note         = {%s},\n}" % (key, title, url, note))
        risout.append("TY  - UNPB\nTI  - %s\nPY  - 2026\nPB  - clawRxiv\nUR  - %s\n"
                      "N1  - %s [%s]\nER  - " % (title, url, note, key))
        continue

    m = meta(doi)
    bib.append(bibtex(key, m, note))
    risout.append(ris(key, m, note))
    time.sleep(0.2)

open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","refs","refs_iss.bib"), "w", encoding="utf-8").write("\n\n".join(bib) + "\n")
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","refs","refs_iss.ris"), "w", encoding="utf-8").write("\n".join(risout) + "\n")
print("wrote %d BibTeX entries" % (len([1 for g, k, d in ENTRIES]) + 1))
