# Part II: Generalizations and Additional Theorems

**Companion to** `Iss_JC69_star_proof.md` (Part I). Part I proved, under JC69 on a star tree: the exact closed form, the second-order form, the $\Gamma$ Laplace transform and its monotonicity in $\alpha$, the non-identifiability of $I_{ss.c}$, and the affine relation to $p$-distance.

Part II removes the JC69/star restrictions, supplies the general criterion for when the affine identity holds, and proves the remaining claims: the two $H_{FSS}$ conventions, the $L$-invariance of the statistic, the Felsenstein/Farris transfer, the nominal range $I_{ss}\le 1$, the uniform-$I_{ss}$ sampling design, and the counting-vector cardinality. Numbering continues from Part I.

---

## S1'. Spectral projector machinery

Throughout, $Q$ is reversible w.r.t. $\boldsymbol\pi$ and normalized by $-\sum_j\pi_jQ_{jj}=1$. Let $\Pi=\mathrm{diag}(\boldsymbol\pi)$ and $\mathbf u_k$ be the orthonormal eigenvectors of the symmetrized generator
$$S=\Pi^{1/2}Q\,\Pi^{-1/2}=U\Lambda U^{\!\top},\qquad \Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_K),\quad 0=\lambda_1>\lambda_2\ge\cdots\ge\lambda_K .$$

**Lemma 9 (spectral resolution).** $\displaystyle P_{ij}(t)=\pi_j+\sum_{k\ge 2}A^{(k)}_{ij}\,e^{\lambda_kt}$, with
$$A^{(k)}_{ij}=\pi_i^{-1/2}u_{ki}\,u_{kj}\,\pi_j^{1/2},\qquad A^{(k)}_{jj}=u_{kj}^2 .$$

*Proof.* $Q=\Pi^{-1/2}U\Lambda U^{\!\top}\Pi^{1/2}$, hence $P(t)=e^{Qt}=\Pi^{-1/2}Ue^{\Lambda t}U^{\!\top}\Pi^{1/2}$ and $P_{ij}(t)=\sum_k\pi_i^{-1/2}u_{ki}e^{\lambda_kt}u_{kj}\pi_j^{1/2}$. For $k=1$: $\lambda_1=0$ and $\mathbf u_1=\sqrt{\boldsymbol\pi}$ (indeed $S\sqrt{\boldsymbol\pi}=\Pi^{1/2}Q\mathbf 1=\mathbf 0$), giving $A^{(1)}_{ij}=\pi_i^{-1/2}\sqrt{\pi_i}\sqrt{\pi_j}\,\pi_j^{1/2}=\pi_j$. $\blacksquare$

**Definition.** For $k\ge2$ set
$$\tau_k:=\operatorname{tr}A^{(k)}=\sum_j u_{kj}^2,\qquad w_k:=\sum_j\pi_j A^{(k)}_{jj}=\sum_j\pi_j u_{kj}^2,\qquad \xi_k:=\big\langle e^{\lambda_k d}\big\rangle\ge 0,$$
where $\langle\cdot\rangle$ averages over ordered leaf pairs (and over the rate distribution when A4 is active). For a distinct eigenvalue $\mu$ define the *grouped* coefficients $\tau_{[\mu]}=\sum_{k:\lambda_k=\mu}\tau_k$ (the multiplicity) and $w_{[\mu]}=\sum_{k:\lambda_k=\mu}w_k$.

**Lemma 10 (the three scalars).** With $D=\sum_j\pi_j\langle P_{jj}\rangle$ and $\wp=1-D$:
$$\Psi=\sum_{k\ge2}\tau_k\,\xi_k,\qquad \wp=\Big(1-\sum_j\pi_j^2\Big)-\sum_{k\ge2}w_k\,\xi_k,\qquad \Delta=\sum_{k\ge2}\Big(w_k-\frac{\tau_k}{K}\Big)\xi_k,$$
where $\Delta:=\sum_j(\pi_j-1/K)\big(\langle P_{jj}\rangle-\pi_j\big)$ is the composition-leakage term of Theorem 4.

*Proof.* By Lemma 9, $\langle P_{jj}\rangle=\pi_j+\sum_{k\ge2}u_{kj}^2\xi_k$. Summing over $j$: $\Psi=\sum_j\langle P_{jj}\rangle-1=1+\sum_{k\ge 2}\tau_k\xi_k-1$. Multiplying by $\pi_j$ and summing: $D=\sum\pi_j^2+\sum_{k\ge2}w_k\xi_k$, and $\wp=1-D$. Finally $\langle P_{jj}\rangle-\pi_j=\sum_{k\ge2}u_{kj}^2\xi_k$, so $\Delta=\sum_j(\pi_j-\tfrac1K)\sum_{k\ge2}u_{kj}^2\xi_k=\sum_{k\ge2}\big(w_k-\tfrac{\tau_k}{K}\big)\xi_k$, using $\sum_j u_{kj}^2=\tau_k$ and $\sum_j\pi_j u_{kj}^2=w_k$. $\blacksquare$

---

## S2'. The two conventions for $H_{FSS}$

**Table S5.** The index is a ratio, so its normalization is a modelling choice with observable consequences.

| | Convention **A** (Xia / DAMBE) | Convention **B** (Part I) |
|---|---|---|
| $H_{FSS}$ | $\mathbb E_{\mathrm{Mult}(N,\hat{\boldsymbol\pi})}[h(\hat{\mathbf n}/N)]$ | $h(\hat{\boldsymbol\pi})$ |
| finite-$N$ bias | **removed** | **retained** |
| leading value | $h(\boldsymbol\pi)-\dfrac{K-1}{2N}+O(N^{-2})$ | $h(\boldsymbol\pi)$ |
| $N=16$, uniform $\pi$ | $1.285617$ nats | $\ln 4=1.386294$ nats |
| ceiling $\lim_{\Psi\to0}\mathbb E[I_{ss}]$ | $1$ | $1-\dfrac{K-1}{2N\,h(\boldsymbol\pi)}=0.932374$ |

*(The $N=16$ entry of column A is the exact multinomial sum $\ln N-\frac1N\sum_j\sum_m m\ln m\binom Nm q_j^m(1-q_j)^{N-m}$ at $q_j=1/4$; the Miller–Basharin term $3/32=0.09375$ nats is only its leading part, the remainder $-0.00693$ nats being $O(N^{-2})$.)*

**Theorem 5 (convention-independence of the affine identity).** Under either convention, $\mathbb E[I_{ss}]$ is an affine function of $\wp$ whenever $\Psi$ is (Theorem 6). The slope and intercept are convention-dependent, but the *criterion* for affinity is not: it depends only on $Q$ and $\boldsymbol\pi$ through the coefficients of Lemma 10.

*Proof.* Both conventions divide $\mathbb E[\bar H]=h(\boldsymbol\pi)-\frac12\big[\frac{K-1}{N}+\frac{N-1}{N}\Psi\big]+O(3)$ by a quantity that does not involve $\Psi$: $h(\hat{\boldsymbol\pi})$ (B) or $\mathbb E_{\mathrm{Mult}(N,\hat{\boldsymbol\pi})}[h]$ (A), both functions of $\hat{\boldsymbol\pi}$ alone. Affinity in $\Psi$ — hence in $\wp$ — is preserved under division by a $\Psi$-free constant. $\blacksquare$

**Corollary S5.1 (numerical reconciliation).** Under JC69 with convention B, $\wp=\tfrac34(1-\xi)$ and $\Psi=3\xi$, so at $N=16$, $d=3.2$ (star, pairwise $=6.4$): $\wp=0.75(1-e^{-8.5333})=0.749853$ and
$$\mathbb E[I_{ss}]=1.442695\cdot\tfrac{15}{16}\cdot 0.749853-0.082021=0.932170,$$
matching the second-order value $0.9322$ and the simulated $0.9273\!\to\!0.9324$ plateau. Here $1.442695=2/\ln 4$ and $0.082021=\tfrac{3}{2\ln4}-1$.

---

## S3'. Master criterion for the affine identity

**Definition (scaling family).** Fix a topology with branch-length vector $\mathbf d^0$, a rate distribution, and scale by $s>0$: $\mathbf d(s)=s\mathbf d^0$. The induced curve is $\boldsymbol\xi(s)\in\mathbb R^{K-1}$. Say the family is *exponentially non-degenerate* if the functions $\{s\mapsto\xi_k(s)\}_{k\ge2}$, taken one representative per distinct eigenvalue, together with the constant function $1$, are linearly independent.

*(This holds generically: for a star tree under $\Gamma$-rates, $\xi_k(s)=(1-2\lambda_k s d^0/\alpha)^{-\alpha}$, and functions with distinct $\lambda_k$ have distinct singularities in the complex plane, hence are linearly independent.)*

**Theorem 6 (criterion).** Assume an exponentially non-degenerate scaling family. Then the following are equivalent:

1. $\Psi$ is an affine function of $\wp$ along the family;
2. $\Delta$ is an affine function of $\Psi$ along the family;
3. there exists a constant $B$ with $\tau_{[\mu]}+B\,w_{[\mu]}=0$ for every distinct non-zero eigenvalue $\mu$.

When (3) holds, $\Psi=B\,\wp-B\big(1-\sum_j\pi_j^2\big)$.

*Proof.* *(1)$\Rightarrow$(3).* Suppose $\Psi=A+B\wp$ for constants $A,B$. By Lemma 10,
$$\sum_k\tau_k\xi_k(s)=A+B\Big(1-\sum\pi_j^2\Big)-B\sum_k w_k\xi_k(s)\ \Longrightarrow\ \sum_k\big(\tau_k+Bw_k\big)\xi_k(s)=\text{const}.$$
Exponential non-degeneracy forces each coefficient to vanish: $\tau_k+Bw_k=0$ for all $k$; summing within a distinct eigenvalue gives (3). *(3)$\Rightarrow$(1)* is the same computation read backwards, giving $\Psi=-B\sum_kw_k\xi_k=B\wp-B(1-\sum\pi^2)$.
*(2)$\Leftrightarrow$(3).* By Lemma 10, $\Delta=\sum_k(w_k-\tau_k/K)\xi_k$ and $\Psi=\sum_k\tau_k\xi_k$. These are affinely related iff $(w_k-\tau_k/K)=c\,\tau_k$ for a constant $c$ (same linear-independence argument), i.e. $w_k=(c+1/K)\tau_k$ for all $k$ — which is (3) with $B=-(c+1/K)^{-1}$. $\blacksquare$

**Remark 10.** Theorem 6 corrects a natural but incorrect reading: **$\Delta\neq0$ does not by itself break the identity.** What breaks it is $\Delta$ failing to be *proportional* to $\Psi$. Both possibilities are realized numerically in Table S7.

**Corollary 6.1 (uniform composition — the strongest case).** If $\boldsymbol\pi$ is uniform, the identity holds for **every** reversible $Q$, **every** topology, **every** branch-length vector, and **every** rate distribution (including $+\Gamma$ and $+I$). Explicitly, with $B=-K$,
$$\Psi=K\wp-\frac{K(K-1)}{K}=3-4\wp\quad(K=4),$$
$$\boxed{\ \mathbb E[I_{ss}]=-\Big(\frac{K-1}{2\ln K}-1\Big)+\frac{2(N-1)}{N\ln K}\,\wp+O(3)\;=\;1.4427\cdot\frac{N-1}{N}\,\wp-0.0820+O(3)\ }$$

*Proof.* Uniform $\pi$ gives $w_k=\frac1K\sum_j u_{kj}^2=\tau_k/K$, so (3) holds with $B=-K$. No condition on the spectrum, topology or rates was used. Substituting into Theorem 1 and collecting constants:
$$1-\frac{1}{2\ln K}\Big[\frac{K-1}{N}-\frac{K(N-1)}{N}\wp+\frac{K(N-1)}{N}\cdot\frac{K-1}{K}\Big]=1-\frac{K-1}{2\ln K}+\frac{K(N-1)}{2N\ln K}\wp,$$
using $\frac{K-1}{N}+\frac{K(N-1)}{N}\cdot\frac{K-1}{K}=\frac{(K-1)K}{N}$. With $K=4$: intercept $1-\frac{3}{2\ln4}=-0.082021$, slope $\frac{2(N-1)}{N\ln4}=1.442695\frac{N-1}{N}$. $\blacksquare$

**Corollary 6.2 (degenerate spectrum).** If $Q$ has a single distinct non-zero eigenvalue $\mu$ (JC69, F81), the identity holds for **every** $\boldsymbol\pi$, with $B=-\tau_{[\mu]}/w_{[\mu]}$.

*Proof.* There is only one grouped pair, so (3) is solvable by $B=-\tau_{[\mu]}/w_{[\mu]}$. For F81, $P(t)=e^{-t}I+(1-e^{-t})\mathbf 1\boldsymbol\pi^{\!\top}$ (since $Q=\mathbf 1\boldsymbol\pi^{\!\top}-I$ and $(\mathbf 1\boldsymbol\pi^{\!\top})^2=\mathbf 1\boldsymbol\pi^{\!\top}$), so $P_{jj}(t)=\pi_j+(1-\pi_j)e^{-t}$, $\mu=-1$, $\tau_{[\mu]}=K-1=3$, $w_{[\mu]}=\sum_j\pi_j(1-\pi_j)=1-\sum\pi_j^2$, and $\wp=(1-\sum\pi_j^2)(1-\langle e^{-d}\rangle)$. $\blacksquare$

**Corollary 6.3 (generic failure).** If $\boldsymbol\pi$ is non-uniform **and** $Q$ has at least two distinct non-zero eigenvalues with $\tau_{[\mu]}/w_{[\mu]}$ not constant across $\mu$, the identity fails: $I_{ss}$ is then a genuinely $m$-dimensional function of $\boldsymbol\xi$ and cannot be recovered from $\wp$. In particular **Xia et al. (2003) simulate under F84 with $\kappa\in[1,10]$ and frequencies in $[0.1,0.9]$** — squarely inside the failure region.

**Lemma 11 (properties of $\Delta$).** $\Delta(0)=\frac1K-\sum_j\pi_j^2$ (all $\xi_k=1$); $\Delta(\infty)=0$; and $\Delta$ need not be monotone, admitting an interior extremum when the coefficients $w_k-\tau_k/K$ change sign.

*Proof.* $\Delta(0)=\sum_k(w_k-\tau_k/K)=\sum_k w_k-\frac{K-1}{K}=\big(1-\sum\pi^2\big)-\frac{K-1}{K}=\frac1K-\sum\pi_j^2$; $\Delta(\infty)=0$ since all $\xi_k\to0$. Non-monotonicity: $\Delta(s)$ is a signed mixture of exponentials with distinct rates. $\blacksquare$

---

## S4'. Rate heterogeneity: general model

**Theorem 7.** Under A4 with $r\sim\Gamma(\alpha,\text{scale}=1/\alpha)$,
$$\Psi_\alpha=\sum_{k\ge2}\tau_k\Big(1-\frac{\lambda_k\bar d}{\alpha}\Big)^{-\alpha},$$
where $\bar d$ enters through $\xi_k=\mathbb E_r[\langle e^{\lambda_k r d}\rangle]$; on a star tree with root-to-leaf length $d$, $\xi_k=(1-2\lambda_k d/\alpha)^{-\alpha}$.

*Proof.* Immediate from Lemma 10 and the Laplace transform $\mathbb E[e^{\lambda rt}]=(1-\lambda t/\alpha)^{-\alpha}$ (Lemma 6), applied term by term; $\lambda_k<0$ for $k\ge2$ guarantees convergence. $\blacksquare$

**Remark 11.** Each spectral term decays as a **power law** $d^{-\alpha}$ rather than exponentially, so *every* mode is retarded. This is the mechanism behind the $\sim\!50$-fold span in the tree length needed to reach a fixed index value (Part I, Remark 5), and it acts on top of — not instead of — the criterion of Theorem 6.

---

## S5'. $L$-invariance of the statistic

**Theorem 8.** Let $\bar H=\frac1L\sum_{i}H_i$ and $B=h(\hat{\boldsymbol\pi})$ with $\hat{\boldsymbol\pi}$ the pooled frequency vector. Then

1. $\mathbb E[\bar H]=\mu(s)$ **exactly**, independent of $L$;
2. $\operatorname{Var}(\bar H)=\sigma^2(s)/L$ with $\sigma^2=\operatorname{Var}(H_1)$;
3. $\operatorname{Var}(I_{ss})=\Theta(L^{-1})$;
4. $\mathbb E[I_{ss}]=I_{ss}^{\infty}+\dfrac{c}{L}+O(L^{-2})$ with $I_{ss}^{\infty}=\mu/h(\boldsymbol\pi)$ and
$$c=\frac{L\,\mu\,\operatorname{Var}(B)}{h(\boldsymbol\pi)^3}-\frac{L\,\operatorname{Cov}(\bar H,B)}{h(\boldsymbol\pi)^2}.$$

*Proof.* (1) Sites are i.i.d. (A5), so $\mathbb E[\bar H]=\mathbb E[H_1]=\mu$. (2) Variance of an i.i.d. mean. (3)–(4) $I_{ss}=\bar H/B$; expanding $B^{-1}$ to second order about $h(\boldsymbol\pi)$,
$$\mathbb E[\bar H B^{-1}]=\frac{\mu}{\mu_B}-\frac{\operatorname{Cov}(\bar H,B)}{\mu_B^2}+\frac{\mu\,\operatorname{Var}(B)}{\mu_B^3}+O(L^{-2}),$$
and $\operatorname{Var}(I_{ss})=\big(\tfrac{\mu}{\mu_B}\big)^2\big[\tfrac{\operatorname{Var}(\bar H)}{\mu^2}-\tfrac{2\operatorname{Cov}(\bar H,B)}{\mu\,\mu_B}+\tfrac{\operatorname{Var}(B)}{\mu_B^2}\big]$, all three terms being $\Theta(L^{-1})$ because $\hat{\boldsymbol\pi}$ is pooled over $NL$ observations. $\blacksquare$

**Corollary 8.1 (the diagnostic contrast).** $I_{ss}$ is $L$-invariant to $O(L^{-1})$ (an intensive, per-site quantity), whereas $I_{ss.c}$ increases monotonically in $L$ (Xia et al. 2003, and §1.2 of the source document). Hence
$$\exists L^*:\quad I_{ss}=I_{ss.c}(L^*),$$
and the verdict flips at $L^*$. Saturation is a property of the *process*; it cannot be reduced by sequencing more sites. That $I_{ss.c}$ tracks $L$ proves it encodes **resolving power**, not saturation.

**Remark 12.** Theorem 8(4) shows the *mean* drifts at $O(L^{-1})$ even though the *median* does not (median is invariant under the monotone per-site structure). Detecting this $c/L$ term empirically is itself independent evidence that the entropy contributes only finite-sample noise.

---

## S6'. The Felsenstein / Farris transfer

Consider the unrooted quartet $((A,B),(C,D))$ with pendant lengths $(a,b,c,d)$ and internal branch $r$.

**Lemma 12.** $\displaystyle\sum_{i<j}d_{ij}=3(a+b+c+d)+4r .$

*Proof.* $d_{AB}=a+b$, $d_{CD}=c+d$, $d_{AC}=a+r+c$, $d_{AD}=a+r+d$, $d_{BC}=b+r+c$, $d_{BD}=b+r+d$. Each pendant index occurs in exactly three of the six sums; $r$ occurs in exactly four. $\blacksquare$

**Definition (the transfer).** Fix $p>q>0$, set $a=p$, $d=q$, and let
$$b=q+t,\qquad c=p-t,\qquad t\in[0,\;p-q].$$
At $t=0$ the configuration is $(p,q,p,q)$: the two long branches $a,c$ are separated by $r$ — the **Felsenstein zone**. At $t=p-q$ it is $(p,p,q,q)$: the long branches $a,b$ are adjacent — the **Farris zone**. The true topology is $((A,B),(C,D))$ throughout.

**Theorem 9.** Along the transfer:

1. **Patristic invariance.** $\sum_{i<j}d_{ij}=3(a+b+c+d)+4r=6p+6q+4r$ is constant, hence so is the mean *tree* distance.
2. **Failure of $p$-distance and $I_{ss}$ invariance.** Mean $p$-distance and $\mathbb E[I_{ss}]$ are *not* constant: they depend on $\sum_{i<j}e^{-4d_{ij}/3}$, which is not a function of $\sum_{i<j}d_{ij}$.
3. **Exact coincidence at $r=0$.** The pairwise-distance multisets of the two endpoints coincide *exactly*:
$$\{\,2p,\;2q,\;(p+q)\times 4\,\},$$
hence $I_{ss}$ is **identical** at the two endpoints.
4. **Gap at order $r$.** For small $r$,
$$\mathbb E[I_{ss}]^{\mathrm{Fels}}-\mathbb E[I_{ss}]^{\mathrm{Farr}}=\frac{N-1}{4N\,h(\boldsymbol\pi)}\cdot\frac{4r}{3}\Big(e^{-4p/3}-e^{-4q/3}\Big)^2+O(r^2)\;>\;0 .$$

*Proof.* (1) Lemma 12 with $a+b+c+d=p+(q+t)+(p-t)+q=2p+2q$.
(3) At $r=0$, Felsenstein gives $\{p+q,\;p+q,\;2p,\;p+q,\;p+q,\;2q\}$ and Farris gives $\{2p,\;2q,\;p+q,\;p+q,\;p+q,\;p+q\}$ — the same multiset. Since $\Psi=\frac12\sum_{i<j}e^{-4d_{ij}/3}$ under uniform $\pi$ (Part I), $\Psi$ — and hence $\mathbb E[I_{ss}]$ — coincide exactly.
(4) Write $S^{\mathrm{F}}(r)=\sum_{i<j}e^{-4d_{ij}/3}$ for the two configurations. Differentiating at $r=0$:
$$\frac{dS^{\mathrm{F}}}{dr}\Big|_0=-\frac43\Big[e^{-8p/3}+e^{-8q/3}+2e^{-4(p+q)/3}\Big],\qquad \frac{dS^{\mathrm{Fa}}}{dr}\Big|_0=-\frac43\Big[4e^{-4(p+q)/3}\Big],$$
so $S^{\mathrm{F}}(r)-S^{\mathrm{Fa}}(r)=-\frac{4r}{3}\big(e^{-4p/3}-e^{-4q/3}\big)^2+O(r^2)$, using $e^{-8p/3}+e^{-8q/3}-2e^{-4(p+q)/3}=(e^{-4p/3}-e^{-4q/3})^2$. Since $\Psi=\frac12 S$ and $\mathbb E[I_{ss}]=1-\frac{1}{2h(\pi)}\big[\frac{K-1}{N}+\frac{N-1}{N}\Psi\big]$, the stated gap follows with $K=4$. $\blacksquare$

**Remark 13 (a correction).** It is tempting to conclude from Lemma 12 that "mean $p$-distance is invariant along the transfer". **This is false**: $p$-distance is $p_{ij}=\frac34(1-e^{-4d_{ij}/3})$, a non-linear function of $d_{ij}$, and Lemma 12 constrains only the *linear* functional. Theorem 9(2)–(4) is the correct statement; the two endpoints agree exactly only in the limit $r\to0$. This matters because $r\to0$ is precisely the regime where the two zones are defined, so the argument is sharpest exactly where it is needed.

**Table S6.** Verification, $N=4$, JC69, $p=1.0$, $q=0.05$. Analytic values from Theorem 1 with $\Psi=\frac12\sum_{i<j}e^{-4d_{ij}/3}$.

| $r$ | Felsenstein | Farris | midpoint | gap | $O(r)$ prediction (Thm 9.4) |
|---|---|---|---|---|---|
| 0.000 | 0.468316 | 0.468316 | 0.501714 | $0.000000$ | $0.000000$ |
| 0.020 | **0.473433** | **0.471826** | **0.505589** | $+0.001607$ | $+0.001628$ |
| 0.050 | 0.480858 | 0.476920 | 0.511211 | $+0.003938$ | $+0.004071$ |

(The $r=0.02$ row reproduces the values 0.4734 / 0.4718 / 0.5056 quoted in the source document to all printed digits.)

**Corollary 9.1 (the experiment).** At fixed $I_{ss}$ (up to $O(r)$), the two configurations have opposite failure modes — long-branch *attraction* versus long-branch *repulsion* — and therefore generically different $\Pr[\hat{\mathcal T}=\mathcal T]$. Plotting $I_{ss}$ against $\Pr[\hat{\mathcal T}=\mathcal T]$ across the transfer (and across reconstruction methods: ML, MP, NJ) yields, at essentially one value of $I_{ss}$, success probabilities spanning the full interval. This simultaneously refutes "saturation causes unrecoverability" and "$I_{ss}$ measures saturation".

**Remark 14.** Because ML is consistent under a correctly specified model, the differences in Corollary 9.1 are finite-$L$ power differences, not inconsistencies — which is exactly the right register, since $I_{ss.c}$ is itself a finite-$L$ construct. To obtain genuine inconsistency one must use MP or NJ with $p$-distances; both should be reported.

---

## S7'. Nominal range: $I_{ss}\le 1$

**Theorem 10.** Under A1–A5 and either convention, $\mathbb E[I_{ss}]\le 1+O(3)+O(N^{-2})$, with equality (convention A) attained only as $\Psi\to0$ (full saturation).

*Proof.* By Lemma 4 and Corollary S1, $\mathbb E[\bar H]=h(\boldsymbol\pi)-\frac12\big[\frac{K-1}{N}+\frac{N-1}{N}\Psi\big]+O(3)$ and $\Psi=\sum_k\tau_k\xi_k\ge0$ (all $\tau_k>0$, all $\xi_k\ge0$). Convention B: divide by $h(\hat{\boldsymbol\pi})=h(\boldsymbol\pi)+O_p((NL)^{-1/2})$. Convention A: $H_{FSS}=\mathbb E_{\mathrm{Mult}(N,\hat{\boldsymbol\pi})}[h]=h(\hat{\boldsymbol\pi})-\frac{K-1}{2N}+O(N^{-2})$, so
$$\mathbb E[I_{ss}]=\frac{h(\boldsymbol\pi)-\frac{K-1}{2N}-\frac{N-1}{2N}\Psi+O(3)}{h(\hat{\boldsymbol\pi})-\frac{K-1}{2N}+O(N^{-2})}\le 1+O(3)+O(N^{-2}),$$
with equality iff $\Psi=0$. $\blacksquare$

**Corollary 10.1.** Observed values $I_{ss}>1$ cannot be produced by any substitution process. They require the numerator and denominator to be evaluated over **different state spaces or different site sets**. The published D-loop case (Iss rising 0.988 → 1.112 → 1.288 → 1.490 as NumOTU goes 4 → 8 → 16 → 32) is explained by: (i) gaps/ambiguity codes counted as a fifth state in the numerator while the denominator stays 4-nucleotide, giving headroom $\ln 5/\ln 4 = 1.160964$; (ii) AT-skew lowering $h(\hat{\boldsymbol\pi})$; (iii) monotonicity in $N$, since $\Pr[\text{site contains a gap}]$ increases with the number of OTUs.

**Remark 15.** Theorem 10 is proved to second order. A fully non-asymptotic version requires comparing $\mathbb E[h]$ under the exchangeable (positively correlated) site law with the multinomial law in convex order; we do not prove this. It is immaterial for the application: observed excursions to $1.49$ exceed any $O(3)$ or $O(N^{-2})$ slack by two orders of magnitude.

---

## S8'. Uniform-$I_{ss}$ sampling design

**Theorem 11.** Let $\mu(s)=\mathbb E[I_{ss}]$ be strictly increasing on $[s_{\min},s_{\max}]$ with $\mu_{\min}=\mu(s_{\min})$, $\mu_{\max}=\mu(s_{\max})$, $Z=\mu_{\max}-\mu_{\min}$. Drawing $s$ from
$$p_s(s)=\frac{\mu'(s)}{Z}\ \mathbf 1_{[s_{\min},s_{\max}]},\qquad F(s)=\frac{\mu(s)-\mu_{\min}}{Z},$$
makes $I_{ss}$ uniform on $[\mu_{\min},\mu_{\max}]$. Sampling is by solving $\mu(s)=\mu_{\min}+UZ$, $U\sim\mathrm{Unif}[0,1]$.

*Proof.* Change of variables: if $Y=\mu(s)$ then $p_Y(y)=p_s(\mu^{-1}(y))/|\mu'(\mu^{-1}(y))|$. Substituting $p_s=\mu'/Z$ gives $p_Y=1/Z$ on $[\mu_{\min},\mu_{\max}]$. The CDF is $F(s)=\int_{s_{\min}}^s\mu'(u)du/Z$, and inversion is the probability-integral transform. $\blacksquare$

**Corollary 11.1.** No pilot run is needed: the normalizing constant is the increment of $\mu$ itself, and $\mu'$ is analytic,
$$\mu'(s)=-\frac{N-1}{2N\,h(\boldsymbol\pi)}\sum_{k\ge2}\tau_k\,\frac{d\xi_k}{ds},$$
with $\frac{d\xi_k}{ds}=\big\langle\lambda_k d^0 e^{\lambda_k s d^0}\big\rangle$ (homogeneous) or $\propto\big(1-\tfrac{\lambda_k s d^0}{\alpha}\big)^{-\alpha-1}$ (under $\Gamma$). Newton converges in a few steps from the analytic derivative.

**Table S7.** Efficiency gain, $N=16$, $\Gamma$-rates, target band $I_{ss}\in[0.30,0.80]$ within a grid spanning $I_{ss}\in[0.10,0.85]$.

| $\alpha$ | $s$-range | span | fraction of uniform-$s$ grid in band | gain |
|---|---|---|---|---|
| 0.25 | $[0.1128,\ 2150.66]$ | $19{,}072\times$ | 14.7% | **6.82×** |
| 0.50 | $[0.0909,\ 28.21]$ | $310\times$ | 37.0% | **2.70×** |
| 1.00 | $[0.0820,\ 4.24]$ | $52\times$ | 54.3% | 1.84× |
| 4.00 | $[0.0760,\ 1.31]$ | $17\times$ | 65.0% | 1.54× |

**Remark 16.** Note the ceiling: under convention B the supremum of $I_{ss}$ is $1-\frac{K-1}{2N\ln K}=0.932374$ ($N=16$), so $I_{ss}=0.95$ is unreachable for *any* $\alpha$. Designs targeting very high index values must either adopt convention A or increase $N$.

---

## S9'. Exact distribution: counting vectors

**Lemma 13.** $H_i$ depends on the site only through the count vector $\mathbf n=(n_1,\dots,n_K)$, $\sum_j n_j=N$. The number of distinct count vectors is $\binom{N+K-1}{K-1}=\binom{N+3}{3}$:
$$\begin{array}{c|ccc} N & 4 & 16 & 32\\\hline \#\mathbf n & 35 & 969 & 6545 \end{array}$$

*Proof.* $h(\hat{\mathbf p})$ is a symmetric function of the leaf states, hence a function of $\mathbf n$ alone. The count of weak compositions of $N$ into $K$ parts is the stars-and-bars number $\binom{N+K-1}{K-1}$. $\blacksquare$

**Proposition 13.1 (polynomial pruning).** On a general tree, attach to each node $v$ the generating polynomial
$$G_v(\mathbf z)=\sum_{\mathbf m}\Pr(\text{counts }\mathbf m\text{ below }v)\,\mathbf z^{\mathbf m},$$
a polynomial in $K$ variables of total degree $\le N_v$. Leaves give $G=z_{X}$; an internal node with children $v_1,v_2$ gives $G_v=\big(P(t_1)\!\cdot\! G_{v_1}\big)\odot\big(P(t_2)\!\cdot\! G_{v_2}\big)$ where the action on the state index is matrix multiplication and $\odot$ denotes convolution in the count index. The root polynomial yields the exact law of $\mathbf n$ in $O(\text{nodes}\cdot K^2\cdot\binom{N+3}{3}^2)$ — seconds for $N=16$.

**Corollary 13.2.** From the exact law one obtains simultaneously: (i) exact $\mu(s)$ and $\sigma^2(s)$, removing the $O(3)$ contamination at low divergence; (ii) exact third moments, hence an Edgeworth correction $-\frac{\gamma_1}{6}(x^2-1)\varphi(x)/\sqrt L$; (iii) a **non-asymptotic** Berry–Esseen bound
$$\sup_x\big|F_L(x)-\Phi(x)\big|\le 0.4748\,\frac{\rho}{\sigma^3\sqrt L},\qquad \rho=\mathbb E|H_1-\mu|^3,$$
valid because $H_1\in[0,\ln K]$ is bounded; (iv) for $N=4$ the exact finite-$L$ law via its 256 atoms and characteristic function $\varphi(t)=\big(\sum_j p_je^{ith_j}\big)^L$, inverted by FFT.

**Remark 17 (three caveats).** (i) At low divergence $H_i$ has a large atom at $0$, so residual plots must be stratified by divergence. (ii) For small $\alpha$ near saturation, required $s$ values are astronomically large (power-law tail) — cap at $0.95\,I_{ss}^{\max}$ and say so. (iii) $\Pr[\hat{\mathcal T}=\mathcal T]$ has **no** closed form (ML tree reconstruction is not analytically solvable); it must be simulated. Stating this explicitly increases, rather than decreases, the credibility of the analytic half.

---

## S10'. Corrections to overstated claims

Nothing in this section weakens the main conclusion; each item replaces a claim that would not survive review with one that does.

**C1 — "Iss.c is not identifiable under $\Gamma$" is too strong.** Table S6 of Part I shows $\alpha=0.5$ reaches $\Pr=0.950$ at $TL=128$: the value exists and is finite, merely displaced. **Corrected:** $I_{ss.c}$ depends strongly and monotonically on $\alpha$ (Lemma 7), while the DAMBE lookup table is calibrated under rate homogeneity; applied to data with $\Gamma$ it is *systematically mis-calibrated*.

**C2 — displacement factors.** From the corrected simulation ($N=4$, $L=1000$, $300$ reps): homogeneous $TL_c\approx7$; $\alpha=4$ crosses $0.95$ in $(8,12)$; $\alpha=1$ in $(16,24)$; $\alpha=0.5$ at $\approx128$. The factors are therefore **$\approx1.4\times\ /\ 2.8\times\ /\ 18\times$**, not $4\times/8\times/>\!30\times$. (The larger figures trace to an earlier evolution-operator implementation that inflated effective tree length by $\approx1.4\times$.)

**C3 — statistical power.** At $\Pr=0.95$ with $300$ replicates, $\mathrm{SE}=\sqrt{0.95\cdot0.05/300}=0.0126$, and $300/300$ successes cannot distinguish $\Pr=1.000$ from $\Pr=0.99$ ($0.99^{300}=0.049$). **Remedy:** raise to $\ge2000$ replicates in the transition region, report "$\Pr>0.99$" rather than "$1.000$", and obtain $TL_c$ by fitting the Gumbel form with a confidence interval rather than by eye.

**C4 — the resampling table does not show a verdict flip.** In the 12S example $I_{ss}\approx0.155$ against a minimum threshold of $0.382$: the verdict never reverses. What the table shows is that the **safety margin drifts by a factor of 2** ($I_{ss.c}^{\mathrm{Asym}}$: $0.768\to0.382$). **Corrected title:** "the decision margin is arbitrary", and either supply a real case with $I_{ss}\in(0.382,0.768)$ or drop the flip claim. Note also that $I_{ss.c}^{\mathrm{Sym}}$ is comparatively stable ($0.800\to0.710$); the collapse is concentrated in the pectinate column, which Xia & Lemey themselves label "extremely unlikely".

**C5 — mean $p$-distance is not invariant along the transfer.** See Remark 13.

**C6 — bidirectional mis-calibration.** Xia's calibration departs from real data in two directions at once: omitting $\Gamma$ (which *delays* the decay of $\Pr$, so the true $TL_c$ is larger and the tabulated $I_{ss.c}$ is too small → under-detection), and using a correctly specified model (real data are mis-specified, which *accelerates* the decay, so the true $TL_c$ is smaller → over-detection). **The sign and magnitude of the net error are not determined**, which is why the test "seems to work" in practice: two biases partially cancel. The defensible formulation is: *an empirically serviceable doubly-mis-calibrated test is lucky, not designed.*

---

## S11'. Numerical validation of Part II

**Table S8 — spectral criterion (Theorem 6).** Ratios $\tau_{[\mu]}/w_{[\mu]}$ must be constant across distinct non-zero eigenvalues $\mu$.

| model | $\boldsymbol\pi$ | # distinct $\mu$ | $\lambda_k$ | $\tau_{[\mu]}$ | $w_{[\mu]}$ | ratios | affine? |
|---|---|---|---|---|---|---|---|
| JC69 | uniform | 1 | $-1.3333^{\times3}$ | $[3]$ | $[0.75]$ | $[4.000]$ | **YES** |
| F81 | $(.4,.3,.2,.1)$ | 1 | $-1.4286^{\times3}$ | $[3]$ | $[0.70]$ | $[4.286]$ | **YES** |
| HKY $\kappa{=}5$ | $(.4,.3,.2,.1)$ | 3 | $-0.6329,-1.6456,-2.1519$ | $[1,1,1]$ | $[.2833,.1500,.2667]$ | $[3.529,6.667,3.750]$ | **NO** |
| HKY $\kappa{=}5$ | uniform | 2 | $-0.5714,-1.7143^{\times2}$ | $[1,2]$ | $[0.25,0.50]$ | $[4.000,4.000]$ | **YES** |
| GTR (arbitrary rates) | uniform | 3 | $-0.6224,-1.2568,-2.1209$ | $[1,1,1]$ | $[.25,.25,.25]$ | $[4,4,4]$ | **YES** |
| GTR (arbitrary rates) | $(.4,.3,.2,.1)$ | 3 | $-0.5418,-1.1829,-2.1807$ | $[1,1,1]$ | $[.1228,.3260,.2512]$ | $[8.147,3.067,3.981]$ | **NO** |

Two rows deserve emphasis. **GTR with uniform $\pi$ and three distinct eigenvalues is affine** — confirming Corollary 6.1 in a case where the spectrum is fully non-degenerate. **F81 with skewed $\pi$ is affine** — confirming Corollary 6.2.

**Table S9 — $\Delta/\Psi$ (Theorem 6, equivalence (2)).** Star tree, leaf–leaf $=2d$, $N=16$.

| $d$ | F81 ($\pi$ skew): $\Delta$ | $\Psi$ | $\Delta/\Psi$ | HKY $\kappa{=}5$ ($\pi$ skew): $\Delta$ | $\Psi$ | $\Delta/\Psi$ |
|---|---|---|---|---|---|---|
| 0.05 | $-0.043344$ | 2.600634 | $-0.016667$ | $-0.040098$ | 2.593328 | $-0.015462$ |
| 0.20 | $-0.028236$ | 1.694154 | $-0.016667$ | $-0.018851$ | 1.716949 | $-0.010980$ |
| 0.50 | $-0.011983$ | 0.718953 | $-0.016667$ | $+0.000349$ | 0.840209 | $+0.000415$ |
| 1.00 | $-0.002872$ | 0.172298 | $-0.016667$ | $+0.005904$ | 0.332736 | $+0.017745$ |
| 2.00 | $-0.000165$ | 0.009896 | $-0.016667$ | $+0.002516$ | 0.081095 | $+0.031019$ |
| 8.00 | $\approx0$ | $\approx0$ | $-0.016667$ | $\approx0$ | $\approx0$ | $+0.033333$ |

For F81 the ratio is **exactly constant** at $w/\tau-1/K=0.7/3-0.25=-0.016667$: $\Delta\neq0$ yet $\Delta\propto\Psi$, so the identity holds (Remark 10). For HKY the ratio sweeps from $-0.0155$ to $+0.0333$ (the limit being $w_{\text{slowest}}/\tau-1/K=0.283333-0.25$), so it fails; $\Delta$ peaks near $d\approx1$, i.e. in the mid-divergence region where $I_{ss.c}$ is calibrated (Lemma 11).

---

## S12'. Summary of Part II

> **S1**–**S3.** Every statistic of interest is an explicit linear functional of the spectral coordinates $\boldsymbol\xi$: $\Psi=\sum\tau_k\xi_k$, $\wp=(1-\sum\pi^2)-\sum w_k\xi_k$, $\Delta=\sum(w_k-\tau_k/K)\xi_k$. **$I_{ss}$ is an affine function of $p$-distance exactly when $\tau_{[\mu]}/w_{[\mu]}$ is constant across distinct eigenvalues** — which holds for *every* reversible model whenever $\pi$ is uniform (Corollary 6.1, verified on GTR with three distinct eigenvalues), and for *every* $\pi$ whenever the spectrum is degenerate (Corollary 6.2, verified on F81). It fails for non-uniform $\pi$ with a non-degenerate spectrum — the HKY/F84+GTR regime Xia simulated in. The correct diagnostic is not $\Delta\neq0$ but $\Delta\not\propto\Psi$.
>
> **S5.** $I_{ss}$ is intensive ($L$-invariant to $O(L^{-1})$) while $I_{ss.c}$ grows with $L$; a critical sample size $L^*$ exists at which the verdict flips. The threshold therefore encodes resolving power.
>
> **S6.** Along the Felsenstein$\to$Farris transfer the mean tree distance is invariant but $I_{ss}$ is not; at $r=0$ the two endpoints' distance multisets coincide *exactly*, so identical $I_{ss}$ accompanies opposite failure modes. The $O(r)$ gap is $\frac{4r}{3}(e^{-4p/3}-e^{-4q/3})^2$, verified numerically.
>
> **S7**–**S8.** $I_{ss}\le1$ in expectation for any substitution process; excursions above $1$ diagnose state-space mismatch, not saturation. And because $\mu(s)$ is analytic and monotone, uniform coverage in $I_{ss}$ is achievable analytically with up to a $6.8\times$ reduction in simulation budget.
>
> **S10.** Six claims are weakened to their defensible forms; none of them is load-bearing for the main conclusion.
