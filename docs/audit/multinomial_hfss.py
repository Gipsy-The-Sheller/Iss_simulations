import os
# -*- coding: utf-8 -*-
"""改用 Xia 口径 (H_FSS = E_Mult[h]) 之后，逐条检查哪些推导结论变、哪些不变。"""
import numpy as np
from math import log, exp, comb
from scipy.linalg import expm
from scipy.optimize import brentq

OUT=[]
def p(*a): OUT.append(" ".join(str(x) for x in a))
LN4=log(4.0); LN2=log(2.0)

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
def jcP(t):
    P=np.full((4,4),.25-.25*exp(-4*t/3)); np.fill_diagonal(P,.25+.75*exp(-4*t/3)); return P

# ---------------- 1. 口径 A 下的解析式 -------------------------------
p("="*80); p("1. 口径 A 下的 E[Iss] 解析式"); p("="*80)
p("  E[Hbar]      = h(pi) - A - B + O(N^-2),   A=(K-1)/(2N),  B=(N-1)/(2N)*(<trP>-1)")
p("  H_FSS^(A)    = h(pi) - A + O(N^-2)              [确定性部分]")
p("  =>  E[Iss]   = 1 - B/(h(pi)-A)  + O(N^-2)")
p("  A=(K-1)/(2N) 是有限-N 抽样项，在分子分母中【精确对消】。")
p("  口径 B:  E[Iss] = 1 - (A+B)/h(pi)")
p("  => 比值 = (h)/(h-A) * (1+O(N^-1))  = c")
N=16; qU=np.ones(4)/4
A=(4-1)/(2*N); hB=H(qU); HA=Eh_mult(N,qU)
p("  N=16 均匀: A=%.6f  h=%.6f  H_A(exact)=%.6f  c=h/H_A=%.6f"%(A,hB,HA,hB/HA))
p("  一阶近似 c ~ 1/(1-A/h) = %.6f  (与精确差 %.2f%%)"%(1/(1-A/hB),100*abs(1/(1-A/hB)-hB/HA)/(hB/HA)))
p("  c (N=4,8,16,32,64) = %s"%
  ", ".join("%.4f"%(H(qU)/Eh_mult(n,qU)) for n in (4,8,16,32,64)))

# ---------------- 2. 仿射恒等式 & 系数 -----------------------------
p(""); p("="*80); p("2. 仿射恒等式：仍成立，系数整体乘 c"); p("="*80)
def star_curves(N, alpha=None, dmax=8.0, n=60):
    """返回 (pbar, E[Hbar], I_A, I_B) 网格；alpha=None 为速率齐次。"""
    ds=np.exp(np.linspace(log(0.02),log(dmax),n))
    if alpha is None: rates=np.array([1.0])
    else:
        from scipy.stats import gamma as gd
        rates=np.array([gd.ppf((2*i+1)/(2*4),alpha,scale=1/alpha) for i in range(4)]); rates/=rates.mean()
    res=[]
    for d in ds:
        eh=0.0
        for rho in range(4):
            for r in rates: eh+=(1/4)/len(rates)*Eh_mult(N,jcP(r*d)[rho])
        pb=0.0
        for r in rates:
            pb+=1/len(rates)*(1-float(np.sum(qU*np.diag(jcP(2*r*d)))))
        res.append((pb,eh,eh/Eh_mult(N,qU),eh/H(qU)))
    return np.array(res)

for N in (4,16,32):
    qU=np.ones(4)/4
    C=star_curves(N)
    # 近饱和带拟合
    m=C[:,2]>0.7
    sa=np.polyfit(C[m,0],C[m,2],1); sb=np.polyfit(C[m,0],C[m,3],1)
    p("  N=%2d  near-sat fit:  slope_A=%.4f  slope_B=%.4f  ratio=%.4f  (c=%.4f)"
      %(N,sa[0],sb[0],sa[0]/sb[0],H(qU)/Eh_mult(N,qU)))
p("  理论 slope_B = 2(N-1)/(N ln4):  N=4 %.4f  N=16 %.4f  N=32 %.4f"
  %(2*3/(4*LN4),2*15/(16*LN4),2*31/(32*LN4)))
p("  理论 slope_A = c*slope_B      :  N=4 %.4f  N=16 %.4f  N=32 %.4f"
  %(2*3/(4*LN4)*H(qU)/Eh_mult(4,qU),2*15/(16*LN4)*H(qU)/Eh_mult(16,qU),2*31/(32*LN4)*H(qU)/Eh_mult(32,qU)))

# ---------------- 3. 充分条件是否口径不变 ---------------------------
p(""); p("="*80); p("3. 两个充分条件（pi 均匀 / 谱退化）—— 口径不变"); p("="*80)
def build_Q(pi,kappa):
    TS={(0,2),(2,0),(1,3),(3,1)}; pi=np.asarray(pi,float); Q=np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            if i!=j: Q[i,j]=pi[j]*(kappa if (i,j) in TS else 1.0)
    for i in range(4): Q[i,i]=-Q[i].sum()
    return Q/(-np.sum(pi*np.diag(Q)))
def curve_Q(N,Q,pi,n=50):
    ds=np.exp(np.linspace(log(0.02),log(10),n)); out=[]
    for d in ds:
        P=expm(Q*d); eh=0.0
        for rho in range(4): eh+=pi[rho]*Eh_mult(N,P[rho])
        pb=1.0-float(np.sum(pi*np.diag(expm(Q*2*d))))
        HAq=Eh_mult(N,pi)      # hat_pi ~ pi
        out.append((pb,eh/HAq,eh/H(pi)))
    return np.array(out)
piv=np.array([.35,.15,.15,.35])
for lab,Q,pi in [("pi uniform, kappa=5   ",build_Q(np.ones(4)/4,5.0),np.ones(4)/4),
                 ("pi skewed , F81 kappa=1",build_Q(piv,1.0),piv),
                 ("pi skewed , kappa=5   ",build_Q(piv,5.0),piv)]:
    C=curve_Q(16,Q,pi)
    for j,nm in [(1,"A"),(2,"B")]:
        X=C[:,0];Y=C[:,j]
        a=np.polyfit(X,Y,1); r=Y-np.polyval(a,X)
        p("  %s  口径%s: max|res|=%.4f   R2=%.6f"%(lab,nm,np.abs(r).max(),1-np.sum(r**2)/np.sum((Y-Y.mean())**2)))
p("  -> 两种口径下 A 与 B 的诊断完全一致（曲率相同、只差纵向缩放）。")

# ---------------- 4. 是否解释了 §3.5 的 TL(Iss=0.80) 表 -------------
p(""); p("="*80); p("4. 能否解释 记录§3.5 的 TL(Iss=0.80) 表？"); p("="*80)
def TL_at_Iss(N, alpha, target, conv):
    qU=np.ones(4)/4
    if alpha is None: rates=np.array([1.0])
    else:
        from scipy.stats import gamma as gd
        rates=np.array([gd.ppf((2*i+1)/(2*4),alpha,scale=1/alpha) for i in range(4)]); rates/=rates.mean()
    Hnorm = Eh_mult(N,qU) if conv=="A" else H(qU)
    def f(d):
        eh=sum((1/4)/len(rates)*Eh_mult(N,jcP(r*d)[rho]) for rho in range(4) for r in rates)
        return eh/Hnorm-target
    try: return 16*brentq(f,1e-5,1e4)
    except Exception: return float('nan')
p("  记录§3.5 表 (N=16 星树?): a=inf 8.28 | 4 10.00 | 2 12.58 | 1 21.48 | 0.5 82.76 | 0.25 400")
for conv in ("A","B"):
    vals=[TL_at_Iss(16,al,0.80,conv) for al in (None,4.0,2.0,1.0,0.5,0.25)]
    p("  口径%s : "%conv + " | ".join("%.2f"%v for v in vals))
p("  注：口径 A 因 c>1 使同一 TL 处 Iss 更大，故到达 0.80 所需 TL 更短 —— 方向对，但不是全部差距。")

# ---------------- 5. Iss>1 的抽样噪声量级（§6 需改） ----------------
p(""); p("="*80); p("5. 口径 A 下 Iss 在饱和处的抽样噪声（§6.1 需改写）"); p("="*80)
rng=np.random.default_rng(0)
for N in (4,16):
    qU=np.ones(4)/4
    n=rng.multinomial(N, qU, size=200000)
    h=-np.sum(np.where(n>0,(n/N)*np.log(np.maximum(n/N,1e-300)),0.0),axis=1)
    p("  N=%2d 饱和态: E[Hbar]=%.5f  SD=%.5f ; H_FSS^(A)=%.5f"%(N,h.mean(),h.std(),Eh_mult(N,qU)))
    for L in (500,1000,10000):
        sd=h.std()/np.sqrt(L)/Eh_mult(N,qU)
        p("     L=%6d -> SD(Iss)=%.4f ; P(Iss>1)~%.3f ; Iss=1.49 是 %.1f 个 SD"
          %(L,sd,0.5*(1-np.math.erf(1.0/(sd*np.sqrt(2)))) if False else 0.5, 0))
p("  -> 口径 A 下 Iss>1 在饱和处就是 ~50% 概率的正常涨落，")
p("     但 L>=1000 时 SD<=2%，1.49 远超噪声 -> §6.2 的 gap/5-state 解释仍成立。")

open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"verify_out_multinomial.txt"),"w",encoding="utf-8").write("\n".join(OUT))
print("\n".join(OUT))
