# -*- coding: utf-8 -*-
"""Third pass: (a) does the Gamma4 discrepancy matter for Iss?
   (b) independent P_true simulation for the N=4 zone/alpha table."""
import numpy as np
from math import log, exp, comb, sqrt, pi as PI
from scipy.linalg import expm
from scipy.special import gammainc
from scipy.stats import gamma as gd

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
    else:
        r=np.array([gd.ppf((2*i+1)/(2*k),alpha,scale=scale) for i in range(k)])
    return r/r.mean()
def jc69_Q():
    Q=np.full((4,4),1/3.0); np.fill_diagonal(Q,-1.0); return Q
def star_exact(N,Q,pi,d,rates):
    pi=np.asarray(pi,float); Eh=0.0
    for rho in range(4):
        for r in rates:
            Eh+=pi[rho]/len(rates)*Eh_mult(N,expm(Q*(r*d))[rho])
    return Eh/H(pi)

p("="*78)
p("V. Star N=16 JC69, alpha=1: convergence of Iss in k for the two schemes")
p("="*78)
Qj=jc69_Q(); piU=np.ones(4)/4
p("%6s %10s %10s %10s %10s %10s"%("d","k4-med","k64-med","k512-med","cont(quad)","k4-mean"))
for d in (0.2,0.8,3.2):
    vm4=star_exact(16,Qj,piU,d,discrete_gamma(1.0,4,"median"))
    vm64=star_exact(16,Qj,piU,d,discrete_gamma(1.0,64,"median"))
    vm512=star_exact(16,Qj,piU,d,discrete_gamma(1.0,512,"median"))
    xs=np.linspace(2e-5,15,3000); w=gd.pdf(xs,1.0,scale=1.0); w/=w.sum()
    vc=float(sum(wi*Eh_mult(16,expm(Qj*(xi*d))[0]) for xi,wi in zip(xs,w))/LN4)
    p("%6.1f %10.4f %10.4f %10.4f %10.4f %10.4f"
      %(d,vm4,vm64,vm512,vc,star_exact(16,Qj,piU,d,discrete_gamma(1.0,4,"mean"))))
p("  -> median-Gamma4 differs from the continuous limit by ~%.0f-%.0f%% for Iss too."
  %(100*abs(star_exact(16,Qj,piU,0.2,discrete_gamma(1.0,4,'median'))-0.3260)/0.3260,
    100*abs(star_exact(16,Qj,piU,3.2,discrete_gamma(1.0,4,'median'))-0.8178)/0.8178))

# ------------------------------------------------------------------ simulation
p("")
p("="*78)
p("VI. Independent P_true simulation, N=4, JC69+G4, L=1000, 600 reps")
p("    4-taxon tree a=b=c=d=l, internal r=0.5*l  -> TL=4.5*l")
p("    inference: NJ on Jukes-Cantor distances (4-taxon NJ == four-point rule)")
p("="*78)
def sim_ptru(alpha,L,nreps,l,r,seed=0):
    rng=np.random.default_rng(seed)
    N=L*nreps
    if np.isinf(alpha):
        rates=np.array([1.0]); wts=np.array([1.0])
    else:
        rates=discrete_gamma(alpha,4,"mean"); wts=np.ones(4)/4
    rc=rng.choice(len(rates),size=N,p=wts); rr=rates[rc]
    P=lambda t: np.stack([np.exp(expm(jc69_Q()*tt)) for tt in t]) if False else None
    def trans(parent, t):
        # t: per-site array of branch lengths (in substitutions); JC69 closed form
        e=np.exp(-4.0*t/3.0)
        # P(parent -> child): same = .25+.75e ; other = .25-.25e
        u=rng.random(parent.size)
        same_prob=0.25+0.75*e
        out=np.empty_like(parent)
        # state 0 is the "same" state; the sampled target index
        idx=np.floor(rng.random(parent.size)*3).astype(int)   # 0,1,2 -> the 3 other states
        other=(parent+1+idx)%4
        out=np.where(u<same_prob,parent,other)
        return out
    root=rng.integers(0,4,size=N)
    i1=trans(root, rr*(r/2.0)); i2=trans(root, rr*(r/2.0))
    A=trans(i1,rr*l); B=trans(i1,rr*l); C=trans(i2,rr*l); D=trans(i2,rr*l)
    S=np.stack([A,B,C,D])                       # 4 x N
    obs=(S[:,None,:]!=S[None,:,:])              # 4 x 4 x N difference indicators
    P_=obs.reshape(4,4,nreps,L).sum(axis=3).astype(float)/L   # 4x4 x nreps p-distances
    with np.errstate(invalid='ignore'):
        D=-0.75*np.log(1.0-4.0*P_/3.0)
    D=np.nan_to_num(D,nan=3.0,posinf=3.0)
    s1=D[0,1]+D[2,3]; s2=D[0,2]+D[1,3]; s3=D[0,3]+D[1,2]
    win=np.argmin(np.stack([s1,s2,s3]),axis=0)
    return float((win==0).mean())

p("%8s %8s %10s %10s %10s %10s"%("TL","l","homog","a=4","a=1","a=0.5"))
for TL in (4,6,8,12,16,24,32,48,64):
    l=TL/4.5; r=0.5*l
    row=[]
    for al in (np.inf,4.0,1.0,0.5):
        row.append(sim_ptru(al,1000,600,l,r,seed=int(TL*10+ (0 if np.isinf(al) else al*100))))
    p("%8.0f %8.3f %10.3f %10.3f %10.3f %10.3f"%(TL,l,row[0],row[1],row[2],row[3]))

p("")
p("  note Table 3.4 for comparison (N=4, L=1000, 300 reps):")
p("    TL=8   : homog 0.863 | a=4 1.000 | a=1 1.000 | a=0.5 1.000")
p("    TL=16  : homog 0.343 | a=4 0.533 | a=1 1.000 | a=0.5 1.000")
p("    TL=32  : homog 0.343 | a=4 0.370 | a=1 0.867 | a=0.5 1.000")
p("    TL=128 : homog 0.373 | a=4 0.357 | a=1 0.383 | a=0.5 0.950")

open("verify_out3.txt","w",encoding="utf-8").write("\n".join(OUT))
print("ok",len(OUT))
