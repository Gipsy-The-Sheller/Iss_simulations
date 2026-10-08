import os
# -*- coding: utf-8 -*-
"""
Independent numerical verification of the analytic claims in
  Iss方法学批判_技术要点汇总..md
  Iss饱和指数_方法学批判与复现记录.md
Results are written to verify_out.txt (UTF-8) to avoid Windows GBK console issues.
"""
import numpy as np
from math import log, exp, comb, sqrt, pi as PI
from scipy.linalg import expm
from scipy.special import gammainc
from scipy.stats import gamma as gamma_dist

OUT = []
def p(*a):
    OUT.append(" ".join(str(x) for x in a))

LN4 = log(4.0)
EULER = 0.5772156649015329

# ---------------------------------------------------------------- entropy
def H(pv):
    pv = np.asarray(pv, float)
    nz = pv > 1e-300
    return float(-np.sum(pv[nz] * np.log(pv[nz])))

# ---------------------------------------------------------------- Q matrices
STATES = "ACGT"
TS = {(0,2),(2,0),(1,3),(3,1)}   # transitions A<->G, C<->T

def build_Q(pi, kappa):
    """HKY-style: q_ij = pi_j * kappa if transition else pi_j ; normalised to unit rate.
       kappa = 1  ->  F81 (q_ij = pi_j), the maximally degenerate spectrum."""
    pi = np.asarray(pi, float)
    Q = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            if i != j:
                Q[i,j] = pi[j] * (kappa if (i,j) in TS else 1.0)
    for i in range(4):
        Q[i,i] = -Q[i].sum()
    mu = -float(np.sum(pi*np.diag(Q)))
    return Q/mu

def jc69_Q():
    Q = np.full((4,4), 1.0/3.0)
    np.fill_diagonal(Q, -1.0)
    return Q

# ---------------------------------------------------------------- exact star-tree E[h]
def Eh_multinomial(N, q):
    """E[h] for n ~ Multinomial(N,q); h = -sum p ln p, p = n/N.
       E[h] = ln N - (1/N) sum_j sum_{m=1}^N m ln m C(N,m) q_j^m (1-q_j)^{N-m}"""
    s = 0.0
    for j in range(len(q)):
        qj = q[j]
        if qj <= 0.0:
            continue
        for m in range(1, N+1):
            s += m*log(m)*comb(N,m)*(qj**m)*((1.0-qj)**(N-m))
    return log(N) - s/N

def star_exact(N, Q, pi, d, rates=(1.0,), probs=None, hfss_denom=None):
    """Iss = E[h] / H_FSS for a star tree, root-to-tip distance d, discrete rates."""
    pi = np.asarray(pi, float)
    probs = probs if probs is not None else [1.0/len(rates)]*len(rates)
    denom = hfss_denom if hfss_denom is not None else H(pi)
    Eh = 0.0
    for rho in range(4):
        for r, w in zip(rates, probs):
            P = expm(Q * (r*d))
            Eh += pi[rho]*w*Eh_multinomial(N, P[rho])
    return Eh/denom

def star_pbar(pi, Q, d, rates=(1.0,), probs=None):
    """mean pairwise p-distance, star tree"""
    pi = np.asarray(pi, float)
    probs = probs if probs is not None else [1.0/len(rates)]*len(rates)
    tot = 0.0
    for r, w in zip(rates, probs):
        P = expm(Q * (r*d))
        # P(X_i = X_j) = sum_rho pi_rho sum_k P_{rho k}^2   (== sum_k pi_k P_kk(2d))
        same = float(np.sum(pi * np.sum(P*P, axis=1)))
        tot += w*(1.0 - same)
    return tot

# ---------------------------------------------------------------- discrete Gamma
def discrete_gamma_rates(alpha, k=4):
    """Yang (1994) equal-probability discrete gamma, mean 1."""
    if np.isinf(alpha):
        return np.ones(k)
    scale = 1.0/alpha
    qs = [i/k for i in range(1,k)]
    # category means of the truncated gamma
    def trunc_mean(a, b):
        # E[X 1(a<X<=b)] for Gamma(shape=alpha, scale=1/alpha)
        f = lambda x: gammainc(alpha+1.0, alpha*x)
        return scale*(f(b) - f(a))
    bounds = [0.0] + [gamma_dist.ppf(q, alpha, scale=scale) for q in qs] + [np.inf]
    r = []
    for i in range(k):
        r.append(trunc_mean(bounds[i], bounds[i+1])*k)
    r = np.array(r)
    return r/r.mean()

# =================================================================
p("="*78)
p("A. Star tree, N=16, uniform pi, JC69 : exact closed form vs note Table 2.8")
p("="*78)
Qj = jc69_Q(); piU = np.ones(4)/4.0
p("%6s %6s %10s %10s %10s %10s" % ("d","alpha","MC(note)","exact(note)","2nd-ord","1st/affine"))
rows = [(0.2,np.inf,0.3968,0.3962),(0.4,np.inf,0.6133,0.6138),
        (0.8,np.inf,0.8143,0.8147),(1.6,np.inf,0.9132,0.9128),
        (3.2,np.inf,0.9273,0.9272),
        (0.2,1.0,0.3462,0.3450),(0.8,1.0,0.6583,0.6579),(3.2,1.0,0.8580,0.8581)]
for d,al,mc,ex in rows:
    if np.isinf(al):
        rates = (1.0,)
    else:
        rates = tuple(discrete_gamma_rates(al,4))
    ex_mine = star_exact(16, Qj, piU, d, rates)
    # 2nd order:  Iss = 1 - 1/(2 ln4) [ 3/N + (N-1)/N * 3*<exp(-4 r d_pair/3)> ]
    x = np.mean([exp(-4.0*(r*2*d)/3.0) for r in rates])   # pair distance = 2d
    second = 1.0 - (1.0/(2*LN4))*(3.0/16 + (15.0/16)*3.0*x)
    # affine (= 2nd order rewritten via pbar)
    pbar = star_pbar(piU, Qj, d, rates)
    aff = 1.44269*(15.0/16)*pbar - 0.08202
    p("%6.1f %6s %10.4f %10.4f %10.4f %10.4f   |mine-exact|=%.4f" %
      (d, al, mc, ex, second, aff, abs(ex_mine-ex)))
    p("       [my exact = %.4f ; my pbar = %.4f]" % (ex_mine, pbar))

p("")
p("--- saturation limits (口径 A vs B), N=16 uniform ---")
hfssA = Eh_multinomial(16, piU)
p("  h(pi) = ln4                     = %.6f  (口径 B denominator)" % LN4)
p("  E_Mult(16,unif)[h]              = %.6f  (口径 A denominator)" % hfssA)
p("  Iss_max A = E[h_sat]/E_Mult[h]  = %.6f" % (Eh_multinomial(16,piU)/hfssA))
p("  Iss_max B = E[h_sat]/ln4        = %.6f   [1-1.0820/N = %.6f]" %
  (Eh_multinomial(16,piU)/LN4, 1.0-1.08202/16))
p("  3/(2*N*ln4) = %.6f ; 3/(2*ln4) = %.6f" % (3.0/(2*16*LN4), 3.0/(2*LN4)))

p("")
p("="*78)
p("B. Gamma Laplace transform table (note Table 2.6): E[e^{-4tr/3}] = (1+4t/(3a))^{-a}")
p("="*78)
p("%6s %6s %12s %12s %12s %12s" % ("alpha","t","analytic","k=4","k=64","k=1024"))
for al,t,an in [(0.5,0.5,0.65465),(1.0,2.0,0.27273),(4.0,2.0,0.12960)]:
    an_mine = (1.0 + 4.0*t/(3.0*al))**(-al)
    vals = []
    for k in (4,64,1024):
        r = discrete_gamma_rates(al,k)
        vals.append(np.mean(np.exp(-4.0*t*r/3.0)))
    p("%6.2f %6.2f %12.5f %12.5f %12.5f %12.5f   (analytic mine %.5f)" %
      (al,t,an,vals[0],vals[1],vals[2],an_mine))
p("  Gamma4 relative error: " + ", ".join(
    "%.1f%%" % (100*abs(np.mean(np.exp(-4*t*discrete_gamma_rates(a,4)/3)) - (1+4*t/(3*a))**(-a))/(1+4*t/(3*a))**(-a))
    for a,t in [(0.5,0.5),(1.0,2.0),(4.0,2.0)]))

p("")
p("="*78)
p("C. Zone experiment numbers (note 4.3): 4 taxa, a=p=1.0 d=q=0.05 r=0.02, N=4 JC69")
p("="*78)
def four_tax_Iss(a,b,c,d,r,N=4):
    dists = [a+b, c+d, a+r+c, a+r+d, b+r+c, b+r+d]
    x = np.mean([exp(-4.0*dd/3.0) for dd in dists])
    pbar = 0.75*(1.0-x)
    Iss = 1.0 - (1.0/(2*LN4))*(3.0/N + (N-1.0)/N*3.0*x)
    return dists, pbar, Iss
pp_, qq_, rr_ = 1.0, 0.05, 0.02
for t,lab in [(0.0,"t=0   (p,q,p,q) Felsenstein"),(0.475,"t=0.475 mid"),
              (0.95,"t=0.95 (p,p,q,q) Farris")]:
    dists, pbar, Iss = four_tax_Iss(pp_, qq_+t, pp_-t, qq_, rr_)
    p("  %-28s sum d_ij=%.4f  pbar=%.5f  Iss=%.4f" % (lab, sum(dists), pbar, Iss))
p("  invariant check: 3(a+b+c+d)+4r = %.4f" % (3*(pp_+(pp_+qq_)+qq_)+4*rr_))
p("  1/d base in affine: 2(N-1)/(N ln4) = %.5f ; intercept 1-3/(2ln4)= %.5f"
  % (2*3.0/(4*LN4), 1-3.0/(2*LN4)))

p("")
p("="*78)
p("D. Affine identity: sufficient conditions (note 2.5)")
p("="*78)
def affinity_test(N, Q, pi, label, ds=None):
    if ds is None:
        ds = np.exp(np.linspace(log(0.03), log(8.0), 40))
    Iss = np.array([star_exact(N,Q,pi,d) for d in ds])
    pb  = np.array([star_pbar(pi,Q,d) for d in ds])
    A = np.vstack([pb, np.ones_like(pb)]).T
    coef, res, *_ = np.linalg.lstsq(A, Iss, rcond=None)
    fit = A@coef
    resid = Iss - fit
    R2 = 1.0 - np.sum(resid**2)/np.sum((Iss-Iss.mean())**2)
    rng = Iss.max()-Iss.min()
    p("  %-46s slope=%8.4f  intercept=%8.4f  R2=%.6f  max|res|=%.5f (%.1f%% of range %.3f)"
      % (label, coef[0], coef[1], R2, np.abs(resid).max(), 100*np.abs(resid).max()/rng, rng))
    return coef

piU2 = np.array([0.25,0.25,0.25,0.25])
piv  = np.array([0.35,0.15,0.15,0.35]); piv = piv/piv.sum()
affinity_test(16, build_Q(piU2, 5.0), piU2, "(i)  pi uniform, HKY kappa=5     [expect affine]")
affinity_test(16, build_Q(piv, 1.0),  piv,  "(ii) pi skewed , F81 (kappa=1)   [expect affine]")
affinity_test(16, build_Q(piv, 5.0),  piv,  "(iii)pi skewed , HKY kappa=5     [expect NOT affine]")
affinity_test(16, build_Q(piv, 2.0),  piv,  "(iv) pi skewed , HKY kappa=2     [expect NOT affine]")
p("  theory slope 2(N-1)/(N ln4) = %.4f, intercept(口径B) = %.4f"
  % (2*15.0/(16*LN4), 1-3.0/(2*LN4)))

p("")
p("="*78)
p("E. Composition-leakage decomposition (note 2.7)")
p("   sum_j u_j - 1  ==  (4 - 4 pbar - 4 sum pi_j^2) - 4*Delta")
p("="*78)
def check_Delta(Q, pi, d):
    P = expm(Q*d); u = np.diag(P)
    pbar = 1.0 - float(np.sum(pi*u))
    lhs = float(np.sum(u)) - 1.0
    aff = 4.0 - 4.0*pbar - 4.0*float(np.sum(pi**2))
    Delta = (aff - lhs)/4.0
    dlt = pi - 0.25
    Delta_direct = float(np.sum(dlt*(u-pi)))
    return lhs, aff, Delta, Delta_direct, pbar
for lab, Q, pi in [("uniform pi, HKY5   ", build_Q(piU2,5.0), piU2),
                   ("skewed  pi, F81    ", build_Q(piv,1.0), piv),
                   ("skewed  pi, HKY5   ", build_Q(piv,5.0), piv)]:
    for d in (0.3, 2.0):
        lhs, aff, D, Dd, pbar = check_Delta(Q, pi, d)
        p("  %s d=%4.1f : sum u-1=%9.6f  affine=%9.6f  Delta=%9.6f (direct %9.6f) diff=%.2e"
          % (lab, d, lhs, aff, D, Dd, abs(D-Dd)))
p("  -> for F81 + skewed pi, Delta != 0 numerically: see row 2 above.")
p("     (the SECOND-order identity still holds, because there Delta is a")
p("      function of pbar; the note's wording 'Delta absorbed' is loose.)")

p("")
p("="*78)
p("F. Jensen bound  Iss <= 1  (note 6.1)")
p("="*78)
p("  口径 B: Iss = Hbar / h(pi_hat) <= 1 pathwise  [Jensen, h concave]  -> rigorous")
p("  口径 A: H_FSS = E_Mult[h] <= h(pi_hat); no ordering vs Hbar -> NOT pathwise")
p("  counterexample under 口径 A: N=4, every site pattern (A,A,C,G)")
pat = np.array([2,1,1,0])/4.0
Hbar = H(pat)
pihat = pat
hfssA_4 = Eh_multinomial(4, pihat)
p("    Hbar=%.6f  h(pi_hat)=%.6f  E_Mult(4,pi_hat)[h]=%.6f" % (Hbar, H(pihat), hfssA_4))
p("    Iss(A) = %.4f  >  1  ->  Iss<=1 is NOT pathwise under 口径 A" % (Hbar/hfssA_4))
p("    Iss(B) = %.4f  <= 1                                " % (Hbar/H(pihat)))

p("")
p("="*78)
p("G. Two 口径: same slope, different intercept?")
p("="*78)
N=16
p("  slope  2(N-1)/(N ln4) = %.5f  (identical in A and B: the (K-1)/N term")
p("                             is constant in <trP>, so it only shifts the intercept)")
p("  intercept B = 1 - 3/(2 ln4)                    = %.5f" % (1-3.0/(2*LN4)))
p("  intercept A = 1 - 3/(2 ln4) + (3/(2N))/ln4     = %.5f"
  % (1-3.0/(2*LN4) + (3.0/(2*N))/LN4))
p("  (i.e. B: Iss_sat=0.9324, A: Iss_sat=1.0000)")

p("")
p("="*78)
p("H. Misc closed forms")
p("="*78)
p("  JC69:  tr P(t)-1 = 3 exp(-4t/3) ; lambda = -4/3 (x3)  [check t=0.5: %.6f vs %.6f]"
  % (np.trace(expm(jc69_Q()*0.5))-1, 3*exp(-4*0.5/3)))
p("  Yang Gamma4 rates, alpha=0.5: %s" % np.round(discrete_gamma_rates(0.5,4),4))
p("  Yang Gamma4 rates, alpha=1.0: %s" % np.round(discrete_gamma_rates(1.0,4),4))
p("  BE constant check: 0.4748 (Shevtsova) -- analytic constant, not re-derived here")
p("  count vectors: C(4+3,3)=%d  C(16+3,3)=%d  C(32+3,3)=%d"
  % (comb(7,3), comb(19,3), comb(35,3)))
p("  ratio of JS/skew note 7.2 -> not checkable without the design script")
p("")
p("  TL_c from note Table 3.4 by linear interpolation on the P_true grid:")
def interp_TL(grid, P, target=0.95):
    for i in range(len(P)-1):
        if P[i] >= target >= P[i+1]:
            return grid[i] + (grid[i+1]-grid[i])*(P[i]-target)/(P[i]-P[i+1])
    return None
g = [6,8,12,16,24,32,48,64,96,128]
cols = {"homog":[1.000,0.863,0.417,0.343,0.380,0.343,0.373,0.397,0.330,0.373],
        "a=4"  :[1.000,1.000,0.793,0.533,0.393,0.370,0.330,0.357,0.307,0.357],
        "a=1"  :[1.000,1.000,1.000,1.000,0.983,0.867,0.563,0.400,0.367,0.383],
        "a=0.5":[1.000,1.000,1.000,1.000,1.000,1.000,0.993,0.997,0.990,0.950]}
base = interp_TL(g, cols["homog"])
for k,v in cols.items():
    t = interp_TL(g, v)
    p("    %-6s TL_c(linear-interp) = %s   ratio vs homogeneous = %s"
      % (k, ("%.2f"%t) if t else "not reached (<128)", ("%.2f"%(t/base)) if t else "-"))
p("    note claims ratios ~1.4 / 2.8 / 18 (via Gumbel fits)")

open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"verify_out.txt"),"w",encoding="utf-8").write("\n".join(OUT))
print("written", len(OUT), "lines")
