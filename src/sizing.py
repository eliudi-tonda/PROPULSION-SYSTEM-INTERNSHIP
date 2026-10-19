"""Component sizing at the SLS take-off design point (uses varcycle.py gas model)."""
import numpy as np
import varcycle as V
R = 287.0
 
def gamma(T, f):
    c = V.cp(T, f)
    return c / (c - R)
 
def flow_area(m, Tt, Pt, M, f=0.0):
    """Area from the compressible flow function, plus an independent rho*V*A check."""
    g = gamma(Tt / (1 + 0.2 * M * M), f)
    for _ in range(5):                       # iterate gamma at the static temperature
        T = Tt / (1 + (g - 1) / 2 * M * M)
        g = gamma(T, f)
    mfp = np.sqrt(g / R) * M * (1 + (g - 1) / 2 * M * M) ** (-(g + 1) / (2 * (g - 1)))
    A = m * np.sqrt(Tt) / (Pt * mfp)         # method 1: flow function
    P = Pt * (T / Tt) ** (g / (g - 1))
    A_check = m / ((P / (R * T)) * M * np.sqrt(g * R * T))   # method 2: rho*V*A
    return A, A_check
 
def size_engine(st, m0, BPR, cool, f, M19, M2=0.60, htr=0.30, V_ref=25.0):
    mc, mb = m0 / (1 + BPR), m0 * BPR / (1 + BPR)
    out = {}
    out['A2'] = flow_area(m0, st['Tt2'], st['Pt2'], M2)
    out['D_fan'] = np.sqrt(4 * out['A2'][0] / (np.pi * (1 - htr ** 2)))
    out['A25'] = flow_area(mc, st['Tt25'], st['Pt25'], 0.50)
    out['A3'] = flow_area(mc, st['Tt3'], st['Pt3'], 0.25)
    rho3 = st['Pt3'] / (R * st['Tt3'])
    out['A_ref'] = (mc / (rho3 * V_ref),) * 2
    out['A4'] = flow_area(mc * (1 - cool + f), st['Tt4'], st['Pt4'], 1.0, st['fb'])
    out['A5'] = flow_area(mc * (1 + f), st['Tt5'], st['Pt5'], 0.45, f)
    out['A8'] = flow_area(mc * (1 + f), st['Tt5'], st['Pt9'], 1.0, f)
    out['A19'] = flow_area(mb, st['Tt13'], st['Pt19'], M19)
    return out