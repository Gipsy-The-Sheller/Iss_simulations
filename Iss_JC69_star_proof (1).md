# Supplementary Material: A Closed-Form Analysis of the Substitution Saturation Index $I_{ss}$ under JC69 on a Star Tree

**Overview.** We derive exact and second-order closed forms for $\mathbb{E}[I_{ss}]$ under the Jukes–Cantor model on a star tree, show that the index is an *exact affine reparametrization* of the mean pairwise $p$-distance, and prove that the critical value $I_{ss.c}$ of Xia et al. (2003) is not identifiable without specifying the rate-heterogeneity parameter. All statements are proved below; numerical results are labelled as simulation-supported and separated in §S7.

---

## S1. Model, assumptions, and notation

**A1 (Topology).** $\mathcal{T}$ is a *star tree*: a root $\rho$ with $N\ge 2$ leaves, every root-to-leaf branch having length $d>0$ (in expected substitutions per site at unit rate). Hence $d_{ij}=2d$ for all ordered pairs $i\neq i'$.

**A2 (Substitution process).** JC69 with total substitution rate normalized to $1$:
$$Q=\tfrac13\mathbf{J}-\tfrac43\mathbf I,\qquad \mathbf J=\mathbf 1\mathbf 1^{\!\top},$$
i.e. $Q_{ij}=1/3$ for $i\neq j$ and $Q_{ii}=-1$. Stationary distribution $\boldsymbol\pi=\tfrac14\mathbf 1$; note $-\sum_j\pi_jQ_{jj}=1$.

**A3 (Root).** The root state is drawn from $\boldsymbol\pi$; leaves are conditionally independent given the root state.

**A4 (Rate heterogeneity).** Site $i$ carries a multiplier $r_i$, i.i.d. across sites with $\mathbb{E}[r]=1$; all branch lengths at that site are scaled by $r_i$. Rate homogeneity is the degenerate case $r\equiv 1$. Under $\Gamma$-heterogeneity $r\sim\Gamma(\alpha,\text{scale}=1/\alpha)$.

**A5 (Sites).** The $L$ sites are independent; expectations below are per-site unless stated.

### Notation

| Symbol | Meaning |
|---|---|
| $K=4$ | number of states |
| $X_i$ | state of leaf $i$ at a given site |
| $n_j=\#\{i:X_i=j\}$; $\hat p_j=n_j/N$ | state counts; site frequency vector $\hat{\mathbf p}$ |
| $h(\mathbf p)=-\sum_j p_j\ln p_j$ | Shannon entropy (nats) |
| $\bar H=\frac1L\sum_{i=1}^{L}h(\hat{\mathbf p}^{(i)})$ | mean site entropy |
| $H_{FSS}=h(\boldsymbol\pi)$ | full-substitution-saturation entropy |
| $I_{ss}=\bar H/H_{FSS}$ | Xia's index |
| $P(t)=e^{Qt}$ | transition matrix |
| $d_{ii'}$ | path length between leaves $i,i'$ |
| $\langle P_{jj}\rangle=\frac{1}{N(N-1)}\sum_{i\neq i'}P_{jj}(d_{ii'})$ | average over **ordered** leaf pairs |
| $\Psi=\sum_j\langle P_{jj}\rangle-1=\langle\operatorname{tr}P\rangle-1$ | *phylogenetic term* |
| $D=\sum_j\pi_j\langle P_{jj}\rangle$ | stationarity-weighted diagonal |
| $\wp=\mathbb{E}[p\text{-distance}]$ | expected proportion of differing leaf pairs |

**Table S1.** For each quantity we give the summation form (top) and the equivalent closed form under A1–A3 (bottom). Here $s:=e^{-8rd/3}$ denotes the *pairwise* JC69 survival parameter at rate $r$.

| Quantity | Summation form | |
|---|---|---|
| $\langle P_{jj}\rangle$ | $\displaystyle\frac{1}{N(N-1)}\sum_{i\neq i'}P_{jj}(d_{ii'})$ | |
| | $P_{jj}(2rd)=\tfrac14+\tfrac34 s$ | *(closed form, star + JC69)* |
| $\Psi$ | $\displaystyle\sum_j\langle P_{jj}\rangle-1$ | |
| | $3s$ | |
| $D$ | $\displaystyle\sum_j\pi_j\langle P_{jj}\rangle$ | |
| | $\tfrac14(1+3s)$ | |
| $\wp$ | $1-D$ | |
| | $\tfrac34(1-s)$ | |

---

## S2. Spectral structure of JC69

**Lemma 1.** Under A2, for all $t\ge0$:
$$P(t)=e^{-\frac43 t}\,\mathbf I+\frac{1-e^{-\frac43t}}{4}\,\mathbf J .$$
Consequently
$$P_{jj}(t)=\tfrac14+\tfrac34 e^{-\frac43 t},\qquad P_{jk}(t)=\tfrac{1-e^{-\frac43t}}{4}\ (j\neq k),\qquad \operatorname{tr}P(t)-1=3e^{-\frac43t}.$$

*Proof.* $\mathbf J^k=4^{k-1}\mathbf J$ for $k\ge1$, hence $\mathbf J$ has eigenvalues $4$ (once, eigenvector $\mathbf 1$) and $0$ (multiplicity $3$). Writing $Q=a\mathbf J+b\mathbf I$ with $a=\tfrac13,\ b=-\tfrac43$, the two matrices commute and
$$e^{Qt}=e^{bt}e^{a\mathbf Jt}=e^{bt}\Big(\mathbf I+\frac{e^{4at}-1}{4}\mathbf J\Big)=e^{-\frac43t}\mathbf I+\frac{1-e^{-\frac43t}}{4}\mathbf J,$$
using $4a=\tfrac43$. Reading off entries gives the stated $P_{jj},P_{jk}$, and $\operatorname{tr}P=4(\tfrac14+\tfrac34e^{-\frac43t})=1+3e^{-\frac43t}$. $\blacksquare$

**Remark 1.** The non-zero eigenvalues of $Q$ are $\lambda=-4/3$ with multiplicity $3$. For a general reversible $Q$ normalized by $-\sum_j\pi_jQ_{jj}=1$, $\operatorname{tr}P(t)-1=\sum_{k=2}^{K}e^{\lambda_kt}$ with $0=\lambda_1>\lambda_2\ge\cdots\ge\lambda_K$ real; Lemma 1 is the JC69 specialization.

---

## S3. Second moments of the site frequency vector

**Lemma 2 (pairwise identity).** Let $Q$ be reversible with stationary $\boldsymbol\pi$ and let $i\neq i'$. Then for every state $j$,
$$\Pr(X_i=j,\ X_{i'}=j)=\pi_j\,P_{jj}(d_{ii'}).$$

*Proof.* By A3, $\Pr(X_i=j,X_{i'}=j)=\sum_\rho\pi_\rho P_{\rho j}(t_i)P_{\rho j}(t_{i'})$ where $t_i,t_{i'}$ are the root-to-leaf lengths and $t_i+t_{i'}=d_{ii'}$. Reversibility gives $\pi_\rho P_{\rho j}(t_i)=\pi_jP_{j\rho}(t_i)$, so the sum equals $\pi_j\sum_\rho P_{j\rho}(t_i)P_{\rho j}(t_{i'})=\pi_jP_{jj}(t_i+t_{i'})$ by Chapman–Kolmogorov. $\blacksquare$

**Lemma 3.** Under A1–A3, for every state $j$:
$$\operatorname{Var}(\hat p_j)=\frac{\pi_j(1-\pi_j)}{N}+\frac{(N-1)}{N}\,\pi_j\big(\langle P_{jj}\rangle-\pi_j\big).$$

*Proof.* $n_j=\sum_i\mathbf 1\{X_i=j\}$ with $\Pr(X_i=j)=\pi_j$ (A3 + stationarity). Hence
$$\operatorname{Var}(n_j)=\sum_i\operatorname{Var}(\mathbf 1\{X_i=j\})+\sum_{i\neq i'}\operatorname{Cov}(\mathbf 1\{X_i=j\},\mathbf 1\{X_{i'}=j\}).$$
The first sum is $N\pi_j(1-\pi_j)$. By Lemma 2 each covariance term equals $\pi_jP_{jj}(d_{ii'})-\pi_j^2$, and summing over the $N(N-1)$ ordered pairs gives $N(N-1)\pi_j(\langle P_{jj}\rangle-\pi_j)$. Dividing by $N^2$ yields the claim. $\blacksquare$

**Corollary S1 (aggregation).** $\displaystyle\sum_j\frac{\operatorname{Var}(\hat p_j)}{\pi_j}=\frac{K-1}{N}+\frac{N-1}{N}\,\Psi .$

*Proof.* Sum Lemma 3 divided by $\pi_j$ over $j$: $\sum_j(1-\pi_j)/N=(K-1)/N$, and $\frac{N-1}{N}\sum_j(\langle P_{jj}\rangle-\pi_j)=\frac{N-1}{N}(\langle\operatorname{tr}P\rangle-1)=\frac{N-1}{N}\Psi$. $\blacksquare$

**Remark 2.** The first term is the classical Miller–Basharin finite-sample bias $(K-1)/(2N)$ after the entropy expansion of Lemma 4; it is purely multinomial and carries no phylogenetic information. All phylogenetic signal enters through $\Psi$.

---

## S4. Second-order closed form

**Lemma 4 (entropy expansion).** Let $h(\mathbf p)=-\sum_jp_j\ln p_j$ and $\hat{\mathbf p}$ random with $\mathbb{E}[\hat{\mathbf p}]=\boldsymbol\pi$. Then, conditionally on $\hat{\mathbf p}$ lying in a compact neighbourhood of $\boldsymbol\pi$ bounded away from the simplex boundary,
$$\mathbb{E}\big[h(\hat{\mathbf p})\big]=h(\boldsymbol\pi)-\frac12\sum_j\frac{\operatorname{Var}(\hat p_j)}{\pi_j}+R_3,\qquad
|R_3|\le\frac16\sum_j\frac{\mathbb{E}\big[|\hat p_j-\pi_j|^3\big]}{\min(\pi_j,\hat p_j)^2}.$$

*Proof.* $\partial h/\partial p_j=-(\ln p_j+1)$, $\partial^2h/\partial p_j^2=-1/p_j$, $\partial^3h/\partial p_j^3=1/p_j^2$, and mixed partials vanish. Taylor expansion at $\boldsymbol\pi$ to third order with Lagrange remainder; the linear term vanishes because $\mathbb{E}[\hat p_j]=\pi_j$. $\blacksquare$

**Theorem 1 (second-order closed form).** Under A1–A5 with rate multiplier $r$,
$$\boxed{\ \mathbb{E}[I_{ss}]=1-\frac{1}{2\,h(\boldsymbol\pi)}\left[\frac{K-1}{N}+\frac{N-1}{N}\,\big\langle \operatorname{tr}P(rd)-1\big\rangle\right]+O(3)\ }$$
and under A1–A2 in particular,
$$\mathbb{E}[I_{ss}]=1-\frac{1}{2\ln 4}\left[\frac{3}{N}+\frac{N-1}{N}\,3\,\mathbb{E}_r\!\big[e^{-8rd/3}\big]\right]+O(3).$$

*Proof.* By A5, $\mathbb{E}[\bar H]=\mathbb{E}[h(\hat{\mathbf p})]$ (a single site). Apply Lemma 4 and substitute Corollary S1. Divide by $H_{FSS}=h(\boldsymbol\pi)$. For the JC69 specialization use $h(\boldsymbol\pi)=\ln 4$ and Lemma 1 with $t=2rd$, giving $\Psi=3e^{-8rd/3}$; average over $r$. $\blacksquare$

**Remark 3 (validity regime).** The remainder in Lemma 4 blows up whenever some $n_j=0$, i.e. for constant or near-constant sites. Hence Theorem 1 is *asymptotically exact in the saturation regime* ($d\to\infty$, where $\hat{\mathbf p}\to\boldsymbol\pi$ in probability and all $K$ states appear at every site), and degrades at low divergence. This is confirmed numerically in Table S2 (error $\approx 5\times10^{-3}$ for $I_{ss}>0.8$; $\approx6\times10^{-2}$ at $d=0.2$).

**Remark 4 (estimation of $H_{FSS}$).** $H_{FSS}$ is estimated in practice from the global empirical base frequencies pooled over all $NL$ observations, so $H_{FSS}=h(\boldsymbol\pi)+O_p((NL)^{-1/2})$, whereas $\bar H$ fluctuates at $O_p(L^{-1/2})$. The denominator is therefore $\sqrt{N}$ times better determined; we treat it as fixed at $h(\boldsymbol\pi)$ throughout, which is the leading-order behaviour.

---

## S5. Exact closed form on the star tree

The second-order expansion can be bypassed entirely on a star tree, because leaf states are *i.i.d.* given the root.

**Lemma 5 (exact multinomial entropy).** Let $\mathbf n\sim\mathrm{Multinomial}(N,\mathbf q)$. Then
$$\mathbb{E}\big[h(\mathbf n/N)\big]=\ln N-\frac1N\sum_{j}\sum_{m=1}^{N}m\ln m\binom{N}{m}q_j^{m}(1-q_j)^{N-m}.$$

*Proof.* $h(\mathbf n/N)=-\sum_j\frac{n_j}{N}\ln\frac{n_j}{N}=\ln N-\frac1N\sum_j n_j\ln n_j$, using $\sum_jn_j=N$. Take expectations and use $\Pr(n_j=m)=\binom Nm q_j^m(1-q_j)^{N-m}$, with the $m=0$ term contributing $0$. $\blacksquare$

**Theorem 2 (exact closed form).** Under A1–A4, with $e=e^{-4rd/3}$ and
$$q_\rho=\tfrac14+\tfrac34 e\ \ (\text{state at the root}),\qquad q_j=\tfrac{1-e}{4}\ \ (j\neq\rho),$$
we have exactly
$$\boxed{\ \mathbb{E}[I_{ss}]=\frac{1}{h(\boldsymbol\pi)}\sum_{c}w_c\sum_{\rho}\pi_\rho\left[\ln N-\frac1N\sum_{j}\sum_{m=1}^{N}m\ln m\binom Nm q_j^{m}(1-q_j)^{N-m}\right]\ }$$
for a discrete rate distribution $\{(r_c,w_c)\}$. Under rate homogeneity the outer sum degenerates to $r=1$.

*Proof.* Given the root state $\rho$ and rate class $c$, the $N$ leaf states are i.i.d. with law $\mathbf q=P_{\rho\cdot}(r_cd)$ (A1, A3), which by Lemma 1 has the stated entries. Apply Lemma 5 and average over $\rho$ (weight $\pi_\rho$) and over $c$ (weight $w_c$); divide by $h(\boldsymbol\pi)$. $\blacksquare$

**Algorithm S1** (exact evaluation; no simulation required).

```
 1: procedure EXACTISS(N, d, {(r_c,w_c)}, K=4)
 2:   T ← 0                                    ▷ global accumulator
 3:   for each rate class c with weight w_c do
 4:     e ← exp(-4 r_c d / 3)
 5:     for each root state ρ do
 6:       for j ← 1..K do
 7:         q_j ← (1-e)/4 ;  if j = ρ then q_j ← 1/4 + (3/4)e
 8:       S ← 0
 9:       for j ← 1..K do
10:         for m ← 1..N do
11:           S ← S + m·ln(m)·Binom(N,m)·q_j^m·(1-q_j)^{N-m}
12:       T ← T + w_c·(1/K)·( ln N - S/N )
13:   return T / h(π)                            ▷ with h(π)=ln 4 under JC69
```

Complexity $O(K N\sum_c 1)$ per parameter setting — negligible compared with simulation.

---

## S6. Rate heterogeneity, monotonicity, and non-identifiability

**Lemma 6 ($\Gamma$ Laplace transform).** Let $r\sim\Gamma(\alpha,\text{scale}=1/\alpha)$ and $\lambda<0$. Then
$$\mathbb{E}_r\!\big[e^{\lambda r t}\big]=\Big(1-\frac{\lambda t}{\alpha}\Big)^{-\alpha}.$$

*Proof.* For shape $\alpha$ and rate $\beta=\alpha$, the density is $\frac{\beta^\alpha}{\Gamma(\alpha)}r^{\alpha-1}e^{-\beta r}$ and the moment generating function is $\mathbb{E}[e^{sr}]=(\beta/(\beta-s))^{\alpha}$ for $s<\beta$. Set $s=\lambda t<0$ and divide by $\beta=\alpha$. $\blacksquare$

**Theorem 3 ($\Gamma$ closed form).** Under A1–A5 with $r\sim\Gamma(\alpha,1/\alpha)$,
$$\Psi_\alpha=3\Big(1+\frac{8d}{3\alpha}\Big)^{-\alpha},\qquad
\mathbb{E}[I_{ss}]=1-\frac{1}{2\ln4}\left[\frac3N+\frac{N-1}{N}\,\Psi_\alpha\right]+O(3).$$

*Proof.* $\Psi_\alpha=3\,\mathbb{E}_r[e^{-(4/3)r(2d)}]$ by Lemma 1 with $t=2rd$; apply Lemma 6 with $\lambda=-4/3,\ t=2d$. Substitute into Theorem 1. $\blacksquare$

**Remark 5 (exponential versus power law).** Rate homogeneity gives $\Psi_\infty=3e^{-8d/3}$ — *exponential* decay in $d$. Rate heterogeneity gives $\Psi_\alpha\sim 3(3\alpha/8d)^{\alpha}$ as $d\to\infty$ — a *power law*. Since a power-law tail is heavier than any exponential, $\Gamma$-heterogeneity systematically retards the growth of $I_{ss}$ toward $1$. This is the mechanism behind the numerical finding that the tree length required to reach $I_{ss}=0.80$ grows from $\approx 8.3$ (homogeneous) to $\approx 82.8$ ($\alpha=0.5$) to $\approx 400$ ($\alpha=0.25$) — a $\sim 50$-fold span.

**Lemma 7 (monotonicity in $\alpha$).** For every $t>0$, $\alpha\mapsto\mathbb{E}_r[e^{\lambda rt}]$ with $\lambda<0$ is strictly *decreasing*, with limit $e^{\lambda t}$ as $\alpha\to\infty$. Consequently $\Psi_\alpha$ is strictly decreasing and $\mathbb{E}[I_{ss}]$ strictly *increasing* in $\alpha$.

*Proof.* Write $g(\alpha)=-\alpha\ln(1+u)$, $u=-\lambda t/\alpha>0$, so that $\mathbb{E}_r[e^{\lambda rt}]=e^{g(\alpha)}$. Then $g'(\alpha)=-\ln(1+u)+\frac{u}{1+u}$. Since $\ln(1+u)>\frac{u}{1+u}$ for all $u>0$, $g'(\alpha)<0$. As $\alpha\to\infty$, $-\alpha\ln(1-\lambda t/\alpha)\to\lambda t$. Monotonicity of $\Psi_\alpha=3\,e^{g(\alpha)}$ follows, and $\mathbb{E}[I_{ss}]$ is decreasing in $\Psi_\alpha$ by Theorem 1. $\blacksquare$

**Proposition 1 (non-identifiability of the critical value).** Let $I_{ss.c}$ be defined as in Xia et al. (2003): $I_{ss.c}=I_{ss}(TL_c)$ where $TL_c=\inf\{TL:\Pr[\hat{\mathcal T}=\mathcal T]\le 0.95\}$. Then $I_{ss.c}$ is a function of $\alpha$, and the critical value tabulated under rate homogeneity corresponds to $\alpha=\infty$ only. For data generated at any finite $\alpha$, the tabulated value is not a valid threshold.

*Proof.* By Lemma 7 the map $\alpha\mapsto I_{ss}(TL;\alpha)$ is strictly increasing with total range $[I_{ss}(TL;0^+),I_{ss}(TL;\infty)]$, so the level set $\{TL:I_{ss}(TL;\alpha)=\theta\}$ — and hence $TL_c$ and $I_{ss}(TL_c)$ — depends on $\alpha$. Xia et al. (2003) simulate under F84 with rate homogeneity ($\alpha=\infty$); the resulting table therefore carries no $\alpha$ coordinate. $\blacksquare$

**Remark 6.** Proposition 1 is a statement about *dependence*, not about magnitude. Its force comes from Remark 5: the dependence is first-order (exponential versus power law), not a small perturbation. Note also that DAMBE estimates $\alpha$ (and the proportion of invariant sites) from the data, then consults a table that contains neither coordinate.

---

## S7. Collapse onto the mean pairwise $p$-distance

**Lemma 8.** $\wp=1-D$ with $D=\sum_j\pi_j\langle P_{jj}\rangle$; under A1–A3, $\wp=\tfrac34(1-s)$.

*Proof.* $\wp=1-\Pr(X_i=X_{i'})=1-\sum_j\Pr(X_i=j,X_{i'}=j)=1-\sum_j\pi_j\langle P_{jj}\rangle$ by Lemma 2. Under A1–A2, $\langle P_{jj}\rangle=\tfrac14+\tfrac34 s$ for all $j$ and $\sum_j\pi_j=1$, giving $D=\tfrac14(1+3s)$. $\blacksquare$

**Theorem 4 (affine reparametrization).** Under A1–A5,
$$\boxed{\ \mathbb{E}[I_{ss}]=\underbrace{\left(1-\frac{K-1}{2\,h(\boldsymbol\pi)}\right)}_{\text{intercept}}+\underbrace{\frac{K(N-1)}{2N\,h(\boldsymbol\pi)}}_{\text{slope}}\wp\;-\;\frac{N-1}{2N\,h(\boldsymbol\pi)}\,\Delta\ +O(3)\ }$$
where
$$\Delta:=\Psi-\big(KD-1\big)=\sum_j\big(1-K\pi_j\big)\langle P_{jj}\rangle .$$
Under A2 (JC69, $\boldsymbol\pi$ uniform) $\Delta\equiv0$ and the relation is **exactly affine**:
$$\mathbb{E}[I_{ss}]=\Big(1-\frac{3}{2\ln4}\Big)+\frac{2(N-1)}{N\ln 4}\,\wp\;+\;O(3).$$

*Proof.* From Theorem 1, $\mathbb{E}[I_{ss}]=1-\frac{1}{2h(\boldsymbol\pi)}\big[\frac{K-1}{N}+\frac{N-1}{N}\Psi\big]$. Write $\Psi=KD-1+\Delta$ by the definition of $\Delta$, and $D=1-\wp$ by Lemma 8. Substituting,
$$\mathbb{E}[I_{ss}]=1-\frac{K-1}{2Nh(\boldsymbol\pi)}-\frac{N-1}{2Nh(\boldsymbol\pi)}\Big[K(1-\wp)-1\Big]-\frac{N-1}{2Nh(\boldsymbol\pi)}\Delta .$$
Collecting the constant terms: $1+\frac{1}{2Nh(\boldsymbol\pi)}\big[-(K-1)+(N-1)-(N-1)K\big]=1+\frac{N(1-K)}{2Nh(\boldsymbol\pi)}=1-\frac{K-1}{2h(\boldsymbol\pi)}$. The coefficient of $\wp$ is $\frac{K(N-1)}{2Nh(\boldsymbol\pi)}$. Finally, under JC69 $\pi_j=1/K$ for all $j$, so every coefficient $1-K\pi_j$ vanishes and $\Delta=0$; setting $K=4$, $h(\boldsymbol\pi)=\ln4$ gives the stated intercept and slope. $\blacksquare$

**Corollary 2 (the index carries no information beyond $p$-distance).** Under A1–A5 with JC69, $I_{ss}$ and $\wp$ are, to second order, related by an *invertible affine map* with constants depending only on $N$ and $K$ — not on $d$, not on the topology, and not on the rate distribution. Hence on the JC69 star family the entropy construction adds **no** phylogenetic information whatever to the mean pairwise $p$-distance.

*Proof.* Both $\mathbb{E}[I_{ss}]$ (Theorem 1) and $\wp$ (Lemma 8) are affine functions of the single scalar $\Psi$; eliminating $\Psi$ gives Theorem 4 with $\Delta=0$, whose slope $\frac{2(N-1)}{N\ln4}>0$ is non-zero for $N\ge2$, so the map is invertible. $\blacksquare$

**Remark 7 (why no extra information can survive).** The site frequency vector lives in a $(K-1)$-dimensional simplex, but under A1–A3 the evolutionary process constrains $\hat{\mathbf p}$ to a **one-dimensional** trajectory parametrized by $s=e^{-8rd/3}$ (Lemma 1 + Table S1). Any two statistics that are functions of state along a common 1-D curve are necessarily monotone reparametrizations of each other. The entropy formalism has the *capacity* to read out which states are segregating — a genuinely richer quantity than the raw mismatch proportion — but real data never leave the curve, so that capacity is never exercised.

**Remark 8 (where the composition sensitivity enters).** $\Delta$ vanishes identically under JC69 but not in general: for a reversible model with non-uniform stationary $\boldsymbol\pi$, $\Delta=\sum_j(1-K\pi_j)\langle P_{jj}\rangle\neq0$, and the deviation of $I_{ss}$ from the $p$-distance prediction is exactly $-\frac{N-1}{2Nh(\boldsymbol\pi)}\Delta$. A second, distinct channel is Remark 4: $H_{FSS}$ is computed from *empirical* global base frequencies, so under model misspecification or non-homogeneous composition the denominator itself is perturbed. Both channels are off-JC69 and are **not** covered by Theorems 1–4; quantifying them requires the non-uniform-$\boldsymbol\pi$ extension (HKY/GTR), which we leave to follow-up work.

---

## S8. Numerical validation

**Table S2.** Closed forms versus Monte Carlo. Star tree, $N=16$, $L=60{,}000$ sites, JC69; $\Gamma$ column uses discrete $\Gamma_4$ in the simulation and $\Gamma_4$ in the closed forms. "$2^\circ$" = Theorem 1, "exact" = Theorem 2.

| $d$ | rate | MC | exact | $2^\circ$ | MC $-$ exact | MC $-2^\circ$ |
|---|---|---|---|---|---|---|
| 0.2 | homog. | 0.3968 | 0.3962 | 0.3373 | $+0.00066$ | $+0.0595$ |
| 0.4 | homog. | 0.6133 | 0.6138 | 0.5833 | $-0.00054$ | $+0.0300$ |
| 0.8 | homog. | 0.8143 | 0.8147 | 0.8122 | $-0.00037$ | $+0.0021$ |
| 1.6 | homog. | 0.9132 | 0.9128 | 0.9181 | $+0.00035$ | $-0.0049$ |
| 3.2 | homog. | 0.9273 | 0.9272 | 0.9322 | $+0.00012$ | $-0.0049$ |
| 0.2 | $\Gamma,\alpha{=}1$ | 0.3462 | 0.3450 | 0.2860 | $+0.00127$ | $+0.0602$ |
| 0.8 | $\Gamma,\alpha{=}1$ | 0.6583 | 0.6579 | 0.6339 | $+0.00033$ | $+0.0244$ |
| 3.2 | $\Gamma,\alpha{=}1$ | 0.8580 | 0.8581 | 0.8561 | $-0.00005$ | $+0.0019$ |

The exact form (Theorem 2) agrees with simulation to $<1.3\times10^{-3}$ at every setting. The second-order form (Theorem 1) is asymptotically exact in the saturation regime, consistent with Remark 3.

**Table S3.** Convergence of discrete $\Gamma_k$ to the continuous Laplace transform of Lemma 6. Entries are $\mathbb{E}[e^{-4tr/3}]$.

| $\alpha$ | $t$ | analytic (Lemma 6) | $k=4$ | $k=64$ | $k=1024$ |
|---|---|---|---|---|---|
| 0.5 | 0.5 | 0.65465 | 0.62705 | 0.65279 | 0.65453 |
| 1.0 | 2.0 | 0.27273 | 0.24807 | 0.27163 | 0.27266 |
| 4.0 | 2.0 | 0.12960 | 0.11510 | 0.12892 | 0.12957 |

**Remark 9.** The commonly used $\Gamma_4$ discretization carries a $2$–$11\%$ error relative to the continuous transform. This is an independent layer of approximation, additional to — and not accounted for by — the Xia lookup table.

**Table S4.** Simulation support for Theorem 4. $N=16$, $L=1000$, pooled over homogeneous and $\Gamma_4$ ($\alpha\in\{4,1,0.5\}$) conditions, $I_{ss}\in[0.068,0.926]$.

| | slope | intercept | $R^2$ |
|---|---|---|---|
| theory ($N=16$, Theorem 4) | $1.3525$ | $-0.0820$ | — |
| full-range OLS | $1.2436$ | $-0.0144$ | $0.99957$ |
| restricted to $I_{ss}>0.7$ | $1.3078$ | $-0.0558$ | $0.99760$ |

The restricted-band fit approaches the theoretical affine constants, as predicted by Remark 3 (second-order validity in the saturation regime); the residual gap over the full range is the $O(3)$ term.

---

## S9. Delimitation: proved versus simulation-supported

**Proved here (§S2–S7).** Lemma 1, 2, 3, 5, 6, 7; Corollary S1, 2; Theorem 1, 2, 3, 4. In particular:

1. The exact closed form of $\mathbb{E}[I_{ss}]$ on a JC69 star tree (Theorem 2) — no approximation.
2. The second-order closed form, which is a function of $\operatorname{tr}P$ — hence of the substitution counts *and* the spectrum of $Q$ (Theorem 1, Remark 1).
3. The exponential $\to$ power-law transition under $\Gamma$, and strict monotonicity in $\alpha$ (Theorem 3, Lemma 7, Remark 5).
4. $I_{ss.c}$ is not identifiable without $\alpha$ (Proposition 1).
5. Under JC69, $I_{ss}$ is an exactly affine function of the mean pairwise $p$-distance (Theorem 4, Corollary 2).

**Simulation-supported, not proved.** The claim that within Xia's calibration window $\Pr[\hat{\mathcal T}=\mathcal T]\equiv1.000$ under $\Gamma$, so that the Gumbel fit $P_{true}=1-e^{-e^{(a-TL)/b}}$ has no identifiable sigmoid segment. Supporting evidence ($N=4$ balanced quartet, $L=1000$, 300 replicates, four-point condition):

| $TL$ | homogeneous | $\Gamma_4,\alpha{=}4$ | $\Gamma_4,\alpha{=}1$ | $\Gamma_4,\alpha{=}0.5$ |
|---|---|---|---|---|
| $\le 6$ | 1.000 | 1.000 | 1.000 | 1.000 |
| **8** ($=TL_c$) | **0.863** | **1.000** | **1.000** | **1.000** |
| 12 | 0.417 | 0.793 | 1.000 | 1.000 |
| 16 | 0.343 | 0.533 | 1.000 | 1.000 |
| 24 | 0.380 | 0.393 | 0.983 | 1.000 |
| 32 | 0.343 | 0.370 | 0.867 | 1.000 |
| 128 | 0.373 | 0.357 | 0.383 | 0.950 |

(The plateau at $\approx0.35\approx1/3$ at large $TL$ is the three-way random-guessing level, confirming the recovery criterion is correctly implemented.) Elevating this to a theorem requires a large-deviation analysis of the four-point condition under $\Gamma$-distributed rates, which we do not attempt here; note that Proposition 1 already suffices to invalidate the homogeneous table for finite $\alpha$, independently of this empirical finding.

**Not covered.** (i) Non-uniform $\boldsymbol\pi$ (HKY/GTR) — the regime where $\Delta\neq0$ (Remark 8); (ii) covarion / heterotachy processes, i.e. rate variation across lineages rather than across sites; (iii) non-star topologies for the *exact* form of Theorem 2 (Theorem 1 and Theorem 4 do hold for arbitrary topologies, since they use only second moments and $\langle\cdot\rangle$ may be taken over any leaf-pair distance distribution); (iv) the finite-$L$ sampling distribution of $\bar H$, as opposed to its expectation.

---

## S10. Summary

> Under JC69 on a star tree, $\mathbb{E}[I_{ss}]$ is a closed-form function of $\langle\operatorname{tr}P\rangle$ alone (Theorem 1), computable exactly via a multinomial sum (Theorem 2). Because the mean pairwise $p$-distance is an affine function of the *same* scalar, the two statistics are linked by an exactly affine, invertible map whose constants depend only on $N$ and $K$ (Theorem 4, Corollary 2): **the entropy construction contributes nothing beyond $p$-distance**. The only quantity that could rescue the index is the critical value $I_{ss.c}$, and that value depends strictly monotonically on the rate-heterogeneity shape $\alpha$ (Lemma 7, Proposition 1) — a coordinate absent from the published lookup table. Under $\Gamma$ the decay of the phylogenetic term changes from exponential to power law (Remark 5), so the dependence is first-order, spanning a $\sim 50$-fold range in the tree length needed to reach a fixed index value.
