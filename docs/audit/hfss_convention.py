import os
# -*- coding: utf-8 -*-
"""H_FSS 口径核对：Xia et al. 2003 Eq.(3) 究竟定义的是哪一个？

Xia Eq.(3):
  H_FSS = sum_{NA,NC,NG,NT} [N!/(NA!NC!NG!NT!)] * PA^NA PC^NC PG^NG PT^NT
          * ( - sum_j p_j log2 p_j ),  p_j = N_j/N,  N = sum N_j
即 MULTINOMIAL 期望熵（口径 A），底数为 2。
两份 md 的 proof 用的是 h(pi)（口径 B）。
"""
import numpy as np
from math import log, exp, comb, log2
OUT = []
def p(*a): OUT.append(" ".join(str(x) for x in a))

LN4 = log(4.0)
def H(v):
    v = np.asarray(v, float); nz = v > 1e-300
    return float(-np.sum(v[nz]*np.log(v[nz])))

def Hfss_A(N, q, base=2):
    """Xia Eq.(3) 的字面实现：multinomial 期望熵，底数 base。"""
    s = 0.0
    for j in range(4):
        qj = q[j]
        if qj <= 0: continue
        for m in range(1, N+1):
            if base == 2:
                term = m*log2(m) if m > 0 else 0.0
            else:
                term = m*log(m)
            s += term*comb(N, m)*(qj**m)*((1-qj)**(N-m))
    return (log2(N) if base == 2 else log(N)) - s/N

p("="*80)
p("1. Xia Eq.(3) 的字面求值  vs  h(pi_hat)   [bits]")
p("="*80)
p("%5s %-28s %10s %10s %9s %9s" % ("N","composition","h(pihat)","H_FSS(A)","gap","gap%"))
for N in (4, 8, 16, 32, 64):
    for lab, q in [("uniform 1/4", np.ones(4)/4),
                   ("AT-rich .35,.15,.15,.35", np.array([.35,.15,.15,.35])),
                   ("strong .7,.1,.1,.1", np.array([.7,.1,.1,.1]))]:
        hb = H(q)/log(2); a = Hfss_A(N, q, 2)
        p("%5d %-28s %10.5f %10.5f %9.5f %8.2f%%"
          % (N, lab, hb, a, hb-a, 100*(hb-a)/hb))
p("  渐近 gap = (K-1)/(2N ln2) bits :")
for N in (4,8,16,32,64):
    p("    N=%3d  (K-1)/(2N ln2) = %.5f bits   [Xia 报的 %.5f]"
      % (N, 3/(2*N*log(2)), N*0+ (H(np.ones(4)/4)/log(2) - Hfss_A(N, np.ones(4)/4, 2))))

p("")
p("="*80)
p("2. 两种口径给同一份数据两个成比例的 Iss（逐样本精确）")
p("   I_A = Hbar/H_FSS(A) , I_B = Hbar/h(pihat)  =>  I_A = c * I_B,  c = h/H_FSS(A)")
p("="*80)
p("%5s %-28s %10s %10s %12s" % ("N","composition","c=h/H_A","Iss_sat(B)","Iss_sat(A)"))
def Eh_mult(N, q):
    s = 0.0
    for qj in q:
        if qj <= 0: continue
        for m in range(1, N+1):
            s += m*log(m)*comb(N,m)*(qj**m)*((1-qj)**(N-m))
    return log(N) - s/N
for N in (4,8,16,32):
    for lab, q in [("uniform", np.ones(4)/4), ("AT-rich", np.array([.35,.15,.15,.35]))]:
        c = H(q)/Hfss_A(N, q, 2) if False else H(q)/Eh_mult(N,q)   # c in nats (Iss is dimensionless)
        # 口径 B 的饱和值 = E[h]/h(pi)
        satB = Eh_mult(N, q)/H(q)
        satA = Eh_mult(N, q)/Eh_mult(N, q)
        p("%5d %-28s %10.5f %10.5f %12.5f" % (N, lab, c, satB, satA))
p("  -> 口径 A 的饱和上限恒为 1.000（定义使然）；口径 B 远低于 1，且 N 越小越低。")

p("")
p("="*80)
p("3. 星树 JC69 均匀：两种口径下的仿射关系（N=16）")
p("   用精确闭式 E[Hbar]，除以各自的分母")
p("="*80)
from scipy.linalg import expm
def jcP(t):
    P = np.full((4,4), .25-.25*exp(-4*t/3)); np.fill_diagonal(P, .25+.75*exp(-4*t/3)); return P
def star_Eh(N, d):
    return sum((1/4)*Eh_mult(N, jcP(d)[rho]) for rho in range(4))
def star_pbar(d):
    return 1.0 - float(np.sum((np.ones(4)/4)*np.diag(jcP(2*d))))
N = 16; qU = np.ones(4)/4
HA = Eh_mult(N, qU); hB = H(qU)
p("  h(pi)=%.6f  H_FSS(A)=%.6f  c=%.6f" % (hB, HA, hB/HA))
p("%8s %10s %10s %10s" % ("d","pbar","E[Iss](A)","E[Iss](B)"))
ds = [0.2,0.4,0.8,1.6,3.2,6.0]
xs, ya, yb = [], [], []
for d in ds:
    eh = star_Eh(N, d); pb = star_pbar(d)
    xs.append(pb); ya.append(eh/HA); yb.append(eh/hB)
    p("%8.1f %10.5f %10.5f %10.5f" % (d, pb, eh/HA, eh/hB))
pa = np.polyfit(xs, ya, 1); pb_ = np.polyfit(xs, yb, 1)
p("  口径 A 拟合: slope=%.4f  intercept=%.4f" % (pa[0], pa[1]))
p("  口径 B 拟合: slope=%.4f  intercept=%.4f" % (pb_[0], pb_[1]))
p("  比值 slopeA/slopeB = %.4f   (应等于 c = %.4f 附近的修正值)" % (pa[0]/pb_[0], hB/HA))
p("  闭式预测 slope_B = 2(N-1)/(N ln4) = %.4f" % (2*15/(16*LN4)))
p("  闭式预测 intercept_B = 1 - 3/(2 ln4) = %.4f" % (1-3/(2*LN4)))
p("  闭式预测 slope_A = slope_B * c = %.4f" % (2*15/(16*LN4)*hB/HA))

p("")
p("="*80)
p("4. 两份 md 的 self-consistency 检验：它们的 MC 用的是哪个口径？")
p("="*80)
p("  汇总 §2.8 / 记录 §4.5 星树 N=16 d=3.2 饱和值：MC = 0.9273, 精确闭式 = 0.9272")
p("  口径 B 的饱和上限（精确）= E[h_sat]/h(pi) = %.4f" % (Eh_mult(16,qU)/H(qU)))
p("  口径 A 的饱和上限 = 1.0000")
p("  -> 二者相差 %.4f；MC 命中 口径 B，且与 口径 A 相差 7.3%%。"
  % (1 - Eh_mult(16,qU)/H(qU)))
p("  另：§5.3(记录) 近饱和带实测斜率 1.3078，理论 口径B = %.4f（差 3%%）、"
  % (2*15/(16*LN4)))
p("      理论 口径A = %.4f（差 10%%）。→ 再证其代码用 口径 B。" % (2*15/(16*LN4)*hB/HA))

p("")
p("="*80)
p("5. DAMBE 打印的 Hmax 反推（汇总 §1.3 的旁证）")
p("="*80)
p("  口径 A 的 Hmax（bits）= h_bits(pihat) - (K-1)/(2N ln2) + O(N^-2)")
for N in (4,8,16,32):
    for hb in (1.95, 1.80):
        p("    N=%2d, h(pihat)=%.2f bits -> Hmax(A) = %.3f bits" % (N, hb, hb - 3/(2*N*log(2))))
p("  汇总称 DAMBE 打印 1.65–1.85 bits（N<=32）-> 与 N=16/32 行的口径 A 相符；")
p("  但 N=4 时口径 A 会低到 ~1.4 bits，故该旁证只能在 N 较大时成立。")

open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"verify_out_hfss.txt"), "w", encoding="utf-8").write("\n".join(OUT))
print("\n".join(OUT))
