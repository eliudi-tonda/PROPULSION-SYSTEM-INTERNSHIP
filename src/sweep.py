"""Figure 1: T4 calibration sweep at SLS take-off (ISA, M0).
Needs cycle.py and massflow.py in the same folder as this file."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)            # works from any working directory / IDE

import numpy as np
import matplotlib
matplotlib.use("Agg")               # saves a PNG; no window is opened
import matplotlib.pyplot as plt
from cycle import cycle
from mass_flow import m0

flight = (0.0, 0.0)                 # SLS, M0
F_TARGET = 116988.0                 # published rating, 26,300 lbf
refined = dict(OPR=32.6, BPR=5.1, FPR=1.60, BoosterPR=1.70, cool_frac=0.12,
               eta_f=0.88, eta_b_lpc=0.86, eta_c=0.87, eta_tH=0.90, eta_tL=0.91,
               eta_burn=0.99, eta_m=0.99, pi_d=0.995, pi_b=0.95, pi_n=0.98)

T4s = np.arange(1500, 1700, 10)
thrust, tsfc = [], []
print(f"{'T4 (K)':>7} {'Thrust (N)':>11} {'TSFC':>8}")
for t4 in T4s:
    r = cycle(flight, T4=float(t4), **refined)
    thrust.append(r['Fsp'] * m0)
    tsfc.append(float(r['TSFC_lb']))
    print(f"{t4:7.0f} {thrust[-1]:11.0f} {tsfc[-1]:8.4f}")

NAVY, RED = "#1F3A5F", "#C0472B"
fig, ax = plt.subplots(figsize=(7.2, 4.0), dpi=200)
ax2 = ax.twinx()
ax.plot(T4s, thrust, color=NAVY, marker="o", ms=4, lw=1.8)
ax2.plot(T4s, tsfc, color=RED, marker="s", ms=4, lw=1.8)
ax.axhline(F_TARGET, color="gray", ls=":", lw=1.3)
ax.axvline(1615, color="gray", ls=":", lw=1.3)
ax.text(1505, 117400, "published target\n116,988 N (26,300 lbf)", fontsize=7.5, color="#555")
ax.set_xlabel("Turbine entry temperature, T4 (K)", fontsize=9.5)
ax.set_ylabel("Thrust (N)", color=NAVY, fontsize=9.5)
ax2.set_ylabel("TSFC (lb/(lbf-hr))", color=RED, fontsize=9.5)
ax.tick_params(colors=NAVY, labelsize=8.5)
ax2.tick_params(colors=RED, labelsize=8.5)
ax.grid(alpha=0.25)
ax.set_title("Calibration sweep: T4 vs. thrust and TSFC at SLS take-off (ISA, M0)\n"
             "Dotted lines mark(T4=1615 K) and the published thrust target",
             fontsize=9.5)
plt.tight_layout()
out = os.path.join(HERE, "calib_sweep.png")
plt.savefig(out, bbox_inches="tight", facecolor="white")
print("saved:", out)