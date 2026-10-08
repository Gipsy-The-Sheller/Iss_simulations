# -*- coding: utf-8 -*-
"""Second pass: chase the Gamma4 discrepancy, test affinity algebraically,
fit Gumbel TL_c, and run a small P_true sanity simulation."""
import numpy as np
from math import log, exp, comb, sqrt
from scipy.linalg import expm
from scipy.special import gammainc
from scipy.stats import gamma as gd
from scipy.optimize import curve_fit, minimize

OUT=[]
def p(*a): OUT.append(" ".join(str(x) for x in a))
LN4=log(4.0)

def H(v):
    v=np.asarray(v,float); nz=v>1e-300
    return float(-np.sum(v[nz]*np.log(v[nz])))

def Eh_mult(N,q):
    s=0.0
    for qj in q:
        if qj<=0: continue
        for m in range(1,N+1):
            s+=m*log(m)*comb(N,m)*(qj**m)*((1-qj)**(N-m))
    return log(N)-s/N

def discrete_gamma(alpha,k=4,mode="mean"):
    if np.isinf(alpha): return np.ones(k)
    scale=1.0/alpha
    if mode=="mean":
        b=[0.0]+[gd.ppf(i/k,alpha,scale=scale) for i in range(1,k)]+[np.inf]
        f=lambda x: gammainc(alpha+1.0,alpha*x)
        r=np.array([scale*(f(b[i+1])-f(b[i]))*k for i in range(k)])
    elif mode=="median":
        r=np.array([gd.ppf((2*i+1)/(2*k),alpha,scale=scale) for i in range(k)])
    elif mode=="midpoint":
        b=[0.0]+[gd.ppf(i/k,alpha,scale=scale) for i in range(1,k)]+[np.inf]
        r=np.array([(b[i]+b[i+1])/2 for i in range(k)])
    return r/r.mean()

p("="*78)
p("I. Which Gamma discretisation did the note use?  target E[e^{-4tr/3}]")
p("="*78)
targets={(0.5,0.5):0.62705,(1.0,2.0):0.24807,(4.0,2.0):0.11510}
p("%6s %6s %10s %12s %12s %12s %12s"%("alpha","t","continuous","k4-mean","k4-median","k4-midpt","note k=4"))
for (a,t),tv in targets.items():
    cont=(1+4*t/(3*a))**(-a)
    v={m:float(np.mean(np.exp(-4*t*discrete_gamma(a,4,m)/3))) for m in ("mean","median","midpoint")}
    v64m=float(np.mean(np.exp(-4*t*discrete_gamma(a,64,"mean")/3)))
    v64med=float(np.mean(np.exp(-4*t*discrete_gamma(a,64,"median")/3)))
    p("%6.2f %6.2f %10.5f %12.5f %12.5f %12.5f %12.5f"%(a,t,cont,v["mean"],v["median"],v["midpoint"],tv))
    p("        k=64 mean %.5f   k=64 median %.5f   (note k=64: see table 2.6)"
      %(v64m,v64med))
p("  note's k=1024 values (0.65453 / 0.27266 / 0.12957) vs continuous (0.65465 / 0.27273 / 0.12960):")
for (a,t) in [(0.5,0.5),(1.0,2.0),(4.0,2.0)]:
    cont=(1+4*t/(3*a))**(-a)
    vm=float(np.mean(np.exp(-4*t*discrete_gamma(a,1024,"mean")/3)))
    vmed=float(np.mean(np.exp(-4*t*discrete_gamma(a,1024,"median")/3)))
    p("    a=%.1f t=%.1f  cont=%.5f  k1024-mean=%.5f  k1024-median=%.5f"%(a,t,cont,vm,vmed))

# --- star-tree exact Iss with the two Gamma discretisations
def build_Q(pi,kappa):
    TS={(0,2),(2,0),(1,3),(3,1)}; pi=np.asarray(pi,float); Q=np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            if i!=j: Q[i,j]=pi[j]*(kappa if (i,j) in TS else 1.0)
    for i in range(4): Q[i,i]=-Q[i].sum()
    return Q/(-np.sum(pi*np.diag(Q)))
def jc69_Q():
    Q=np.full((4,4),1/3.0); np.fill_diagonal(Q,-1.0); return Q
def star_exact(N,Q,pi,d,rates):
    pi=np.asarray(pi,float); Eh=0.0
    for rho in range(4):
        for r in rates:
            Eh+=pi[rho]/len(rates)*Eh_mult(N,expm(Q*(r*d))[rho])
    return Eh/H(pi)

p("")
p("="*78)
p("II. Star N=16, alpha=1: which Gamma4 reproduces the note's 'exact' column?")
p("="*78)
Qj=jc69_Q(); piU=np.ones(4)/4
p("%6s %12s %12s %12s %12s"%("d","note exact","k4-mean","k4-median","continuous"))
for d,ex in [(0.2,0.3450),(0.8,0.6579),(3.2,0.8581)]:
    vm=star_exact(16,Qj,piU,d,discrete_gamma(1.0,4,"mean"))
    vmed=star_exact(16,Qj,piU,d,discrete_gamma(1.0,4,"median"))
    # continuous: quadrature over the gamma density
    xs=np.linspace(1e-4,12,400); w=gd.pdf(xs,1.0,scale=1.0); w/=w.sum()
    vc=sum(wi*star_exact(16,Qj,piU,d,np.array([xi])) for xi,wi in zip(xs,w))
    p("%6.1f %12.4f %12.4f %12.4f %12.4f"%(d,ex,vm,vmed,vc))

p("")
p("="*78)
p("III. Affine identity, tested ALGEBRAICALLY on  X=<trP>-1  vs  pbar")
p("     (avoids the O(N^-2) Taylor noise of the exact E[Iss])")
p("="*78)
def test_affine(Q,pi,label):
    ds=np.exp(np.linspace(log(0.02),log(10),60))
    Xs=[];pbs=[]
    for d in ds:
        P=expm(Q*d); u=np.diag(P)
        Xs.append(float(np.sum(u))-1.0)
        pbs.append(1.0-float(np.sum(pi*u)))
    Xs=np.array(Xs);pbs=np.array(pbs)
    A=np.vstack([pbs,np.ones_like(pbs)]).T
    c,*_=np.linalg.lstsq(A,Xs,rcond=None)
    res=Xs-A@c
    p("  %-42s slope=%8.4f int=%8.4f  max|res|=%.2e  R2=%.8f"
      %(label,c[0],c[1],np.abs(res).max(),1-np.sum(res**2)/np.sum((Xs-Xs.mean())**2)))
    return c
piU2=np.array([.25]*4); piv=np.array([.35,.15,.15,.35])
test_affine(build_Q(piU2,1.0),piU2,"pi uniform, JC69/F81 (kappa=1)")
test_affine(build_Q(piU2,5.0),piU2,"pi uniform, kappa=5       -> expect EXACT line")
test_affine(build_Q(piU2,20.),piU2,"pi uniform, kappa=20      -> expect EXACT line")
test_affine(build_Q(piv,1.0),piv, "pi skewed, kappa=1 (F81)  -> expect EXACT line")
test_affine(build_Q(piv,5.0),piv, "pi skewed, kappa=5        -> expect CURVE")
test_affine(build_Q(piv,20.),piv, "pi skewed, kappa=20       -> expect CURVE")

p("")
p("="*78)
p("IV. Gumbel fit of the note's Table 3.4 -> TL_c and the multipliers")
p("="*78)
g=np.array([1,3,5,6,8,12,16,24,32,48,64,96,128],float)
cols={"homog":np.array([1,1,1,1,.863,.417,.343,.380,.343,.373,.397,.330,.373]),
      "a=4"  :np.array([1,1,1,1,1,.793,.533,.393,.370,.330,.357,.307,.357]),
      "a=1"  :np.array([1,1,1,1,1,1,1,.983,.867,.563,.400,.367,.383]),
      "a=0.5":np.array([1,1,1,1,1,1,1,1,1,.993,.997,.990,.950])}
def gumbel(TL,a,b): return 1.0-np.exp(-np.exp((a-TL)/b))
def TLc_from_gumbel(TL,P):
    best=None
    for a0 in np.linspace(-4,10,40):
        for b0 in np.logspace(log(0.5),log(30),40):
            try:
                pr,_=curve_fit(gumbel,TL,P,p0=[a0,b0],maxfev=4000)
                r=np.sum((gumbel(TL,*pr)-P)**2)
                if best is None or r<best[0]: best=(r,pr)
            except Exception: pass
    return best
res={}
for k,v in cols.items():
    r,pr=TLc_from_gumbel(g,v)
    TLc=pr[0]-pr[1]*log(-log(1-0.95))
    res[k]=TLc
    p("  %-6s gumbel a=%.3f b=%.3f  TLc=%.2f  (SSE=%.4f)  linear-interp=%.2f"
      %(k,pr[0],pr[1],TLc,r,
        (lambda vv: next((g[i]+(g[i+1]-g[i])*(vv[i]-.95)/(vv[i]-vv[i+1]) for i in range(len(vv)-1) if vv[i]>=.95>=vv[i+1]),float('nan')))(v)))
p("  ratios to homogeneous: a=4 %.2fx   a=1 %.2fx   a=0.5 %.2fx"
  %(res["a=4"]/res["homog"],res["a=1"]/res["homog"],res["a=0.5"]/res["homog"]))
p("  note claims 1.4x / 2.8x / 18x")

open("verify_out2.txt","w",encoding="utf-8").write("\n".join(OUT))
print("ok",len(OUT))
