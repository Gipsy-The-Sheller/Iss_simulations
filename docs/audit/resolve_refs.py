# -*- coding: utf-8 -*-
"""Resolve every work named in the two .md files to a DOI via Crossref, then
pull BibTeX + RIS by content negotiation.  Writes refs_raw.json + a review log."""
import os, json, time, urllib.parse, urllib.request, re, sys, os

UA = "IssCritiqueBib/1.0 (mailto:researcher@example.org)"
BASE = "https://api.crossref.org/works"

def get(url, accept=None):
    req = urllib.request.Request(url, headers={"User-Agent": UA,
          **({"Accept": accept} if accept else {})})
    with urllib.request.urlopen(req, timeout=40) as r:
        return r.read().decode("utf-8", "replace")

def search(q, rows=6):
    url = BASE + "?" + urllib.parse.urlencode(
        {"query.bibliographic": q, "rows": rows,
         "select": "DOI,title,author,issued,container-title,type,publisher,volume,page,issue"})
    try:
        d = json.loads(get(url))
    except Exception as e:
        return []
    return d["message"]["items"]

def score(it, entry):
    s = 0.0
    t = (it.get("title") or [""])[0].lower()
    for w in re.findall(r"[a-z]{4,}", entry.get("title_words", "").lower()):
        if w in t: s += 1.0
    au = it.get("author") or []
    fams = [a.get("family", "").lower() for a in au]
    for h in entry.get("authors", []):
        if h.lower() in fams: s += 3.0
    yr = entry.get("year")
    if yr:
        iy = (it.get("issued", {}).get("date-parts") or [[None]])[0][0]
        if iy == yr: s += 2.0
        elif iy and abs(iy - yr) <= 1: s += 1.0
        elif iy: s -= 1.5
    ct = (it.get("container-title") or [""])[0].lower()
    for j in entry.get("journals", []):
        if j.lower() in ct: s += 2.5
    if entry.get("want_type") and it.get("type") != entry["want_type"]: s -= 1.0
    return s

ENTRIES = [
 # ---------------- A. 核心：Xia / DAMBE ----------------
 dict(key="xia2003", known_doi="10.1016/S1055-7903(02)00326-3",
      q="An index of substitution saturation and its application Xia",
      title_words="index substitution saturation application", authors=["xia"],
      year=2003, journals=["molecular phylogenetics and evolution"]),
 dict(key="xia_xie2001", q="DAMBE data analysis in molecular biology and evolution 2001",
      title_words="DAMBE data analysis molecular biology evolution",
      authors=["xia"], year=2001, journals=["journal of heredity"]),
 dict(key="xia2013_dambe5", q="DAMBE5 a comprehensive software package for data analysis in molecular biology and evolution",
      title_words="DAMBE5 comprehensive software package data analysis molecular biology evolution",
      authors=["xia"], year=2013, journals=["molecular biology and evolution"]),
 dict(key="xia2018_dambe7", q="DAMBE7 new and improved tools for data analysis in molecular biology and evolution",
      title_words="DAMBE7 new improved tools data analysis molecular biology evolution",
      authors=["xia"], year=2018, journals=["molecular biology and evolution"]),
 dict(key="xia_starless", q="starless bias parameter estimation bias likelihood phylogenetic method Xia",
      title_words="starless bias parameter estimation phylogenetic", authors=["xia"],
      journals=[]),
 dict(key="xia2020_book", q="A Mathematical Primer of Molecular Phylogenetics Xia",
      title_words="mathematical primer molecular phylogenetics", authors=["xia"],
      year=2020, journals=[]),
 dict(key="kubatko_review_xia2020", q="book review A Mathematical Primer of Molecular Phylogenetics",
      title_words="book review mathematical primer molecular phylogenetics",
      authors=["kubatko"], journals=["systematic biology"]),
 dict(key="phypa2016", q="PhyPA phylogenetic method with pairwise sequence alignment",
      title_words="phypa phylogenetic pairwise sequence alignment", authors=["xia"],
      year=2016, journals=[]),

 # ---------------- B. 替代指标 ----------------
 dict(key="struck2008_Cvalue", q="Struck 2008 saturation 18S 28S Annelida phylogenetic reconstruction",
      title_words="saturation annelida phylogenetic reconstruction", authors=["struck"],
      year=2008, journals=[]),
 dict(key="kuck_struck2014_bacoca", q="BaCoCa a heuristic software tool for the parallel assessment of sequence biases",
      title_words="BaCoCa heuristic software tool parallel assessment sequence biases",
      authors=["kuck","struck"], year=2014, journals=[]),
 dict(key="struck2014_trespex", q="TreSpEx detection of misleading signal in phylogenetic reconstructions based on tree information",
      title_words="TreSpEx misleading signal phylogenetic reconstructions tree information",
      authors=["struck"], year=2014, journals=[]),
 dict(key="jermiin_misof2020", q="Jermiin Misof 2020 historical signal lambda variable sites",
      title_words="historical signal phylogenetics", authors=["jermiin","misof"],
      year=2020, journals=[]),
 dict(key="fleming_struck2023_nrcfv", q="nRCFV dataset-size-independent metric compositional heterogeneity nucleotide amino acid",
      title_words="nRCFV dataset size independent metric compositional heterogeneity",
      authors=["fleming","struck"], year=2023, journals=[]),
 dict(key="struck2024_saturation", q="SatuRation saturation index compositional heterogeneity Struck",
      title_words="saturation", authors=["struck"], journals=[]),
 dict(key="satute2025", q="SatuTe branch-wise substitution saturation test phylogenetic Manuel",
      title_words="saturation test branch", authors=["manuel"], year=2025,
      journals=["molecular biology and evolution"]),
 dict(key="duchene2021_phylomad", q="Duchene 2021 model adequacy phylogenetic analysis Systematic Biology",
      title_words="model adequacy phylogenetic", authors=["duchene"], year=2021,
      journals=["systematic biology"]),
 dict(key="duchene2018_phylomad_soft", q="PhyloMAd efficient assessment of phylogenomic model adequacy",
      title_words="PhyloMAd assessment phylogenomic model adequacy",
      authors=["duchene"], journals=[]),
 dict(key="townsend2007_pi", q="Profiling phylogenetic informativeness Townsend",
      title_words="profiling phylogenetic informativeness", authors=["townsend"],
      year=2007, journals=["systematic biology"]),

 # ---------------- C. 先例 / 综述 ----------------
 dict(key="verbruggen2007", q="Verbruggen 2007 Molecular Phylogenetics and Evolution substitution saturation",
      title_words="", authors=["verbruggen"], year=2007,
      journals=["molecular phylogenetics and evolution"]),
 dict(key="bmcevolbiol2011_rrna", q="substitution saturation is not a site-intrinsic property phylogenetic level rRNA 2011 BMC",
      title_words="substitution saturation phylogenetic", year=2011,
      journals=["bmc evolutionary biology"]),
 dict(key="amb_apriori_review", q="a priori assessment of data quality for phylogenetic analysis algorithms for molecular biology",
      title_words="data quality phylogenetic", journals=["algorithms for molecular biology"]),
 dict(key="robust_phylogenetics2025", q="Robust Phylogenetics circularity tree-dependent methods fast-evolving sites",
      title_words="robust phylogenetics circularity", year=2025, journals=[]),

 # ---------------- D. 软件 ----------------
 dict(key="iqtree2", q="IQ-TREE 2 new models and efficient methods for phylogenetic inference in the genomic era",
      title_words="IQ-TREE 2 new models efficient methods phylogenetic inference genomic era",
      authors=["minh"], year=2020, journals=["molecular biology and evolution"]),
 dict(key="alisim", q="AliSim a fast and versatile phylogenetic sequence simulator for the genomic era",
      title_words="AliSim fast versatile phylogenetic sequence simulator genomic era",
      authors=["ly-trong","minh"], year=2022, journals=["systematic biology"]),
 dict(key="paml4", q="PAML 4 phylogenetic analysis by maximum likelihood Yang",
      title_words="PAML 4 phylogenetic analysis maximum likelihood", authors=["yang"],
      year=2007, journals=["molecular biology and evolution"]),
 dict(key="phylobayes", q="PhyloBayes MPI phylogenetic reconstruction with infinite mixtures of profile Dirichlet processes Lartillot",
      title_words="PhyloBayes MPI phylogenetic reconstruction infinite mixtures",
      authors=["lartillot"], year=2009, journals=[]),
 dict(key="beast1.10", q="Bayesian phylogenetic and phylodynamic data integration using BEAST 1.10",
      title_words="bayesian phylogenetic phylodynamic data integration BEAST",
      authors=["suchard"], year=2018, journals=[]),
 dict(key="treetime", q="TreeTime maximum-likelihood phylodynamic analysis Sagulenko",
      title_words="TreeTime maximum likelihood phylodynamic analysis", authors=["sagulenko"],
      year=2018, journals=[]),
 dict(key="nextstrain", q="Nextstrain real-time tracking of pathogen evolution Hadfield",
      title_words="Nextstrain real time tracking pathogen evolution", authors=["hadfield"],
      year=2018, journals=[]),
 dict(key="trad_xia", q="Xia virus substitution rate dating divergence time estimator TRAD",
      title_words="virus substitution rate dating", authors=["xia"], journals=[]),

 # ---------------- E. 其他被点名的 ----------------
 dict(key="mbe40th", q="Molecular Biology and Evolution 40th anniversary software phylogenetic",
      title_words="anniversary", journals=["molecular biology and evolution"]),
 dict(key="mafft", q="MAFFT multiple sequence alignment software version 7 Katoh Standley",
      title_words="MAFFT multiple sequence alignment software version 7", authors=["katoh"],
      year=2013, journals=[]),
 dict(key="muscle", q="MUSCLE multiple sequence alignment with high accuracy and high throughput Edgar",
      title_words="MUSCLE multiple sequence alignment high accuracy high throughput",
      authors=["edgar"], year=2004, journals=[]),
 dict(key="phyml", q="New algorithms and methods to estimate maximum-likelihood phylogenies assessing the performance of PhyML 3.0",
      title_words="new algorithms methods estimate maximum likelihood phylogenies PhyML",
      authors=["guindon"], year=2010, journals=[]),
]

out = {}
log = []
for e in ENTRIES:
    if e.get("known_doi"):
        best = dict(DOI=e["known_doi"], _score="known")
        cands = []
    else:
        cands = search(e["q"])
        scored = sorted(((score(c, e), c) for c in cands), key=lambda x: -x[0])
        best = scored[0][1] if scored else None
        if best is not None:
            best = dict(best); best["_score"] = round(scored[0][0], 2)
    out[e["key"]] = dict(entry=e, chosen=best,
                         top=[dict(DOI=c.get("DOI"),
                                   t=(c.get("title") or [""])[0][:90],
                                   y=(c.get("issued", {}).get("date-parts") or [[None]])[0][0],
                                   j=(c.get("container-title") or [""])[0][:45],
                                   sc=round(s, 2))
                              for s, c in (sorted(((score(c, e), c) for c in cands),
                                                  key=lambda x: -x[0])[:3])])
    log.append("%-24s %-6s -> %s  %s" % (
        e["key"], best.get("_score") if best else "NONE",
        (best or {}).get("DOI", "-"),
        ((best or {}).get("title") or [""])[0][:70] if best else "*** NOT FOUND ***"))
    time.sleep(0.4)

open("refs_raw.json", "w", encoding="utf-8").write(json.dumps(out, ensure_ascii=False, indent=1))
open("refs_resolve_log.txt", "w", encoding="utf-8").write("\n".join(log))
print("\n".join(log))
