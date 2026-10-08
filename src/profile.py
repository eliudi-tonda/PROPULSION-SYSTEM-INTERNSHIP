"""Plot total temperature and total pressure through the core gas path
at the calibrated SLS take-off point. Values are computed directly from
cycle.py """

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from cycle import cycle
 
flight = (0.0, 0.0)
refined = dict(OPR=32.6, BPR=5.1, FPR=1.60, BoosterPR=1.70, cool_frac=0.12,
               eta_f=0.88, eta_b_lpc=0.86, eta_c=0.87, eta_tH=0.90, eta_tL=0.91,
               eta_burn=0.99, eta_m=0.99, pi_d=0.995, pi_b=0.95, pi_n=0.98)
r = cycle(flight, T4=1615.0, **refined)
 
stations = ['0/2\n(ambient/\nfan face)', '13\n(fan exit)', '25\n(booster\nexit)', '3\n(HPC exit)',
            '4\n(burner\nexit)', '45\n(HPT exit)', '5\n(LPT exit)', '9\n(core\nnozzle)']
Tt = [r['T0'], r['Tt13'], r['Tt25'], r['Tt3'], r['Tt4'], r['Tt45'], r['Tt5'], r['core']['T']]
Pt = [r['P0']/1000, r['Pt13']/1000, r['Pt25']/1000, r['Pt3']/1000,
      r['Pt4']/1000, r['Pt45']/1000, r['Pt5']/1000, r['Pt9']/1000]
 
NAVY="#1F3A5F"; RED="#C0472B"
x = np.arange(len(stations))
fig, ax1 = plt.subplots(figsize=(8.6,4.0), dpi=200)
ax2 = ax1.twinx()
l1,=ax1.plot(x, Tt, color=NAVY, marker='o', lw=2, ms=6, label="Total temperature Tt")
l2,=ax2.plot(x, Pt, color=RED, marker='s', lw=2, ms=6, label="Total pressure Pt")
ax2.set_yscale('log')
ax1.set_xticks(x); ax1.set_xticklabels(stations, fontsize=8)
ax1.set_ylabel("Tt (K)", color=NAVY, fontsize=10)
ax2.set_ylabel("Pt (kPa, log scale)", color=RED, fontsize=10)
ax1.tick_params(colors=NAVY); ax2.tick_params(colors=RED)
ax1.grid(alpha=0.25)
ax1.set_title("Core-stream total temperature and pressure through the engine\nSLS take-off, ISA, M0, T4=1615 K (calibrated)", fontsize=10.5)
fig.legend(handles=[l1,l2], loc='upper center', bbox_to_anchor=(0.5,0.02), ncol=2, fontsize=9, frameon=False)
plt.tight_layout(rect=[0,0.04,1,1])
plt.savefig("profile.png", bbox_inches="tight", facecolor="white")
 
if __name__ == "__main__":
    print("Station data plotted (from cycle.py, not hardcoded):")
    for s, tt, pt in zip(stations, Tt, Pt):
        print(f"  {s.replace(chr(10),' '):30s} Tt={tt:7.1f} K   Pt={pt:8.1f} kPa")
    print("saved profile.png")