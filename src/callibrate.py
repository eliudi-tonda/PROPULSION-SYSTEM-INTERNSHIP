"""Two-parameter calibration (T4 and a common efficiency offset) to published thrust and ICAO fuel flow."""
from scipy.optimize import fsolve
import varcycle as V
 
F_PUB  = 116988.0          # N   (26,300 lbf)
MF_PUB = 1.221             # kg/s, ICAO take-off fuel flow
M0     = 752 * 0.45359237  # kg/s, published airflow (752 lb/s)
BPR, OPR = 5.1, 27.2
BASE = dict(eta_f=0.88, eta_b=0.86, eta_c=0.87, eta_tH=0.90, eta_tL=0.91)
 
def residuals(x):
    T4, d = x
    eta = {k: v + d for k, v in BASE.items()}
    r = V.cycle(T4, OPR=OPR, BPR=BPR, **eta)
    return [r['Fsp'] * M0 - F_PUB,                 # thrust residual, N
            r['f'] * M0 / (1 + BPR) - MF_PUB]      # fuel-flow residual, kg/s
 
T4, d_eta = fsolve(residuals, [1600.0, 0.03])
eta = {k: v + d_eta for k, v in BASE.items()}
r = V.cycle(T4, OPR=OPR, BPR=BPR, **eta)
print(f"T4 = {T4:.1f} K   eta offset = {d_eta:+.4f}   TSFC = {r['TSFC']:.4f} lb/(lbf.h)")
print(f"Tt3 = {r['Tt3']:.1f} K  Tt45 = {r['Tt45']:.1f} K  Tt5 = {r['Tt5']:.1f} K  f = {r['f']:.5f}")
 
if __name__ == "__main__":
    import sizing, numpy as np
    st = {k: r[k] for k in ('Tt2','Pt2','Tt13','Pt13','Tt25','Pt25','Tt3','Pt3','Tt4','Pt4','Tt5','Pt5','Pt9','Pt19','fb')}
    PR = r['Pt19'] / 101325.0; g = sizing.gamma(r['Tt13'] * 0.97, 0.0)
    M19 = np.sqrt(2 / (g - 1) * (PR ** ((g - 1) / g) - 1))
    s = sizing.size_engine(st, M0, BPR, r['cool'], r['f'], M19)
    for k, v in s.items():
        print(k, np.round(v, 4))