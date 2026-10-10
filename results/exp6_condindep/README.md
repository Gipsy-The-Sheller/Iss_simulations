# Experiment 6 — conditional independence

**Question.** After controlling for the mean pairwise p-distance, does $I_{ss}$
still carry information about a target $Y$?

$$\Pr(Y\mid \bar p, I_{ss})\quad\text{vs}\quad\Pr(Y\mid \bar p)$$

**Design.** Crossed, because scale alone is not enough: if the tree is scaled
with everything else fixed, $I_{ss}$ and $\bar p$ are deterministically linked
and the test is vacuous.  Here internal-branch fraction (4) × pendant
arrangement (2: equal / graded 1:3) × composition and model (5) × scale (18)
are crossed, 150 replicates each, L = 1000, JC69-driven AliSim.
**108,000 scored alignments.**

## Result

Three targets, all evaluated with a **cubic spline in $\bar p$** as the
baseline (see Method note below).

| target $Y$ | range | $\bar p$ only | $+\,I_{ss}$ | Δ |
|---|---|---:|---:|---:|
| MPBoot recovery (AUC) | full | 0.8886 | 0.8895 | **+0.0008** |
| ML recovery, IQ-TREE (AUC) | full | 0.9143 | 0.9163 | **+0.0020** |
| parsimony-informative fraction (CV $R^2$) | full | 0.2268 | 0.2966 | **+0.0698** |
| MPBoot recovery (AUC) | band | 0.5552 | 0.7236 | **+0.1684** |
| ML recovery (AUC) | band | 0.6036 | 0.7851 | **+0.1815** |
| parsimony-informative fraction | band | 0.0450 | 0.2551 | **+0.2101** |

*Band* = $\bar p\in[0.63,0.70]$, where recovery is neither 0 nor 1.  In the
band, $I_{ss}$ **alone** (no $\bar p$ term) reaches AUC 0.750 (MPBoot) and
0.800 (ML), and a single threshold at $I_{ss}\approx0.98$ classifies 80% / 83%
correctly.

### Method note: the baseline must be flexible

The original analysis used a degree-2 polynomial in $\bar p$ and reported
Δ = −0.0103 for MPBoot recovery.  That number is an artefact.  In the
transition band a quadratic in $\bar p$ has cross-validated AUC 0.37 — *worse
than chance* — because the fit is non-monotone and inverts the ordering there;
against such a baseline, adding $I_{ss}$ looks harmful over the full range
while it is in fact repairing a mis-specified competitor.  With a cubic spline
the full-range effect for recovery is 0.000, not a loss, and the band effect is
+0.17.

## Reading

**The second-order argument does not exhaust what $I_{ss}$ contains.**  Beyond
second order $\E[I_{ss}]$ depends on third and fourth moments of the
site-pattern distribution, which pairwise identity does not determine.  $I_{ss}$
predicts the parsimony-informative fraction well beyond $\bar p$ (CV $R^2$
0.227 → 0.297 over the full range), so "whatever $I_{ss}$ adds beyond
$p$-distance is composition" is false as a general statement and must be scoped
to second order.

**$I_{ss}$ is not redundant for recovery either — but only in a narrow band at
the ceiling.**  Over the full range $\bar p$ alone already separates almost
perfectly (AUC ≈ 0.89–0.91) and $I_{ss}$ adds nothing.  In the band where the
decision is actually made it adds $+0.17$ to $+0.18$ AUC.  Two things bound
that:

1. **The window is a sliver.**  Across the band — which spans the whole
   transition from reliable to unreliable recovery — $I_{ss}$ varies only
   between 0.79 and 1.03 (IQR 0.87–0.97).  The index discriminates recovery
   only within a few percent of its ceiling, which is exactly the regime
   §2.6 identifies as one where it cannot discriminate *divergence*.  The
   information is real, but it is not saturation information.
2. **Part of it is tree shape.**  The internal-branch fraction — a property of
   the true tree that no user observes — reaches AUC 0.759 in the band for
   MPBoot (better than $I_{ss}$'s 0.724) and raises the full-range AUC from
   0.889 to 0.923, after which $I_{ss}$ adds nothing.  For ML the ordering
   reverses (0.753 vs 0.785), and in all three targets $I_{ss}$ retains
   something after shape is included, so the two are not substitutes.

## Caveats

1. Recovery here is **MPBoot** (Hoang et al. 2018), a published program.  Both
   the four-point rule and an in-process Fitch implementation were dropped: the
   project reports topology inference only from peer-reviewed tools.  The
   four-point run gave Δ = −0.0092 and MPBoot gives Δ = −0.0103, so the
   conclusion is not an artefact of the criterion.  **ML recovery remains
   unrun**: it costs one IQ-TREE call per alignment (~0.7 s), i.e. ~21 h for all
   108,000, or ~35 min on a 3,000-alignment subsample.
2. The informativeness fraction is a *pattern-composition* statistic, so the
   $R^2$ gain partly reflects that both it and $I_{ss}$ are functionals of the
   site-pattern distribution while $\bar p$ is a first-moment statistic.  The
   gain therefore measures *access to pattern structure*, not necessarily
   *useful* extra information.
3. Only 4 taxa; the higher-order channel is expected to shrink with $N$ but this
   was not tested.
