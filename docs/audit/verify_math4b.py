# -*- coding: utf-8 -*-
"""Is the note's steep Gamma plateau (P_true=1.000 up to TL=32) a property of
   the DATA or of the INFERENCE MODEL?  Compare NJ(JC), homogeneous ML and
   Gamma-aware ML on exactly the same simulated alignments."""
import numpy as np
from math import log, exp
from scipy.optimize import minimize

OUT=[]
def p(*a): OUT.append(" ".join(str(x) for x in a))

# ---------- exact 4-taxon likelihood via pattern counts (unrooted ((A,B),(C,D)))
PERM=np.array([[a,b,c,d] for a in range(4) for b in range(4) for c in range(4) for d in range(4)])
def jcP(t):
    e=np.exp(-4.0*np.asarray(t)/3.0)
    P=np.empty(t.shape+(4,4))
    P[...]=0.25-0.25*e[...,None,None]
    i=np.arange(4)
    P[...,i,i]=0.25+0.75*e[...,None]
    return P
def gamma_rates(alpha,k=4):
    if np.isinf(alpha): return np.array([1.0]),np.array([1.0])
    from scipy.stats import gamma as gd
    from scipy.special import gammainc
    sc=1.0/alpha
    b=[0.0]+[gd.ppf(i/k,alpha,scale=sc) for i in range(1,k)]+[np.inf]
    f=lambda x: gammainc(alpha+1.0,alpha*x)
    r=np.array([sc*(f(b[i+1])-f(b[i]))*k for i in range(k)])
    return r/r.mean(), np.ones(k)/k

def loglik_one(counts_mask, br, rates, wts):
    """counts_mask: (256,)+shape of index positions for this topology; br: 5 lengths."""
    tot=0.0
    L=np.zeros(256)
    pi=np.ones(4)/4
    for r,w in zip(rates,wts):
        P=jcP(r*np.array(br))          # (5,4,4)
        M=[P[i] for i in range(5)]     # A,B,internal,C,D  (M[i][x,y]=P(x->y))
        # L(a,b,c,d)= sum_{u,v} pi_u Ma[u,a] Mb[u,b] Mint[u,v] Mc[v,c] Md[v,d]
        T=np.einsum('ua,ub,uv,vc,vd->abcd', M[0],M[1],M[2],M[3],M[4])
        T*=pi[:,None,None,None]
        L+=w*T.ravel()
    L=np.maximum(L,1e-300)
    return float(np.dot(counts_mask,L))

TOPO_PERM=[PERM, PERM[:, [0,2,1,3]], PERM[:, [0,3,1,2]]]

def sim_counts(alpha,L,nreps,l,r,seed):
    rng=np.random.default_rng(seed); N=L*nreps
    rates,wts=gamma_rates(alpha)
    rc=rng.choice(len(rates),size=N,p=wts); rr=rates[rc]
    def tr(par,t):
        e=np.exp(-4.0*t/3.0); u=rng.random(par.size)
        idx=np.floor(rng.random(par.size)*3).astype(int)
        return np.where(u<0.25+0.75*e,par,(par+1+idx)%4)
    root=rng.integers(0,4,size=N)
    i1=tr(root,rr*r/2); i2=tr(root,rr*r/2)
    A=tr(i1,rr*l); B=tr(i1,rr*l); C=tr(i2,rr*l); D=tr(i2,rr*l)
    ids=(A*64+B*16+C*4+D)
    cts=np.zeros((nreps,256),int)
    for i in range(nreps):
        cts[i]=np.bincount(ids[i*L:(i+1)*L],minlength=256)
    return cts, np.stack([A,B,C,D])

def ml_ptru(cts, states, inferred_alpha, nreps):
    rates,wts=gamma_rates(inferred_alpha)
    wins=np.zeros(nreps,int)
    for i in range(nreps):
        ll=[]
        for t,per in enumerate(TOPO_PERM):
            idx=per[:,0]*64+per[:,1]*16+per[:,2]*4+per[:,3]
            cm=cts[i][idx]
            f=lambda x: -loglik_one(cm,np.exp(x),rates,wts)
            best=minimize(f,[log(0.5)]*5,method="Nelder-Mead",
                          options=dict(maxiter=250,xatol=2e-3,fatol=1e-3))
            ll.append(-best.fun)
        wins[i]=int(np.argmax(ll))
    return float((wins==0).mean())

def nj_ptru(states, L, nreps):
    S=states.reshape(4,nreps,L)
    obs=(S[:,None,:,:]!=S[None,:,:,:])
    Pp=obs.sum(axis=3).astype(float)/L
    with np.errstate(invalid='ignore',divide='ignore'):
        D=-0.75*np.log(1.0-4.0*Pp/3.0)
    D=np.nan_to_num(D,nan=3.0,posinf=3.0)
    s=np.stack([D[0,1]+D[2,3],D[0,2]+D[1,3],D[0,3]+D[1,2]])
    return float((np.argmin(s,axis=0)==0).mean())

p("="*78)
p("VII. Same simulated data, three inference criteria (N=4, L=1000, 300 reps)")
p("     tree: pendants l, internal r=0.5*l, TL=4.5*l  ; JC69+G4")
p("="*78)
p("%6s %8s | %9s %9s %9s | %9s %9s"%("TL","alpha","NJ(JC)","ML homog","ML+G4","note",""))
for TL,alpha in [(16,0.5),(32,0.5)]:
    l=TL/4.5; r=0.5*l
    cts,st=sim_counts(alpha,1000,200,l,r,seed=int(TL*100+(0 if np.isinf(alpha) else alpha*10)))
    nj=nj_ptru(st,1000,200)
    mlh=ml_ptru(cts,st,np.inf,200)
    mlg=ml_ptru(cts,st,alpha,200)
    p("%6.0f %8s | %9.3f %9.3f %9.3f"%(TL,alpha,nj,mlh,mlg))
    p("   (note Table 3.4: TL=%d alpha=%s -> see below)"%(TL,alpha))
p("")
p("  note Table 3.4:")
p("    TL=8  : a=inf 0.863 | a=1 1.000 | a=0.5 1.000")
p("    TL=16 : a=inf 0.343 | a=1 1.000 | a=0.5 1.000")
p("    TL=32 : a=inf 0.343 | a=1 0.867 | a=0.5 1.000")
p("    TL=64 : a=inf 0.397 | a=1 0.400 | a=0.5 0.997")
p("    TL=128: a=inf 0.373 | a=1 0.383 | a=0.5 0.950")
open("verify_out4.txt","w",encoding="utf-8").write("\n".join(OUT))
print("ok")
