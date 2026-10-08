"""Preliminary SLS static mass flow estimate from the fan annulus,
used only to convert specific thrust into a Newtons figure for
comparison"""

import numpy as np

R, gam = 287.0, 1.4
D_tip = 1.549      # m, CFM56-7B26 published fan diameter
htr = 0.30         # assumed hub-to-tip ratio, typical first pass
M_face = 0.60      # assumed fan-face axial Mach number, typical design value

A = (np.pi/4)*D_tip**2*(1-htr**2)

T0, P0 = 288.15, 101325.0   # SLS ISA static = total at M0=0
rho0 = P0/(R*T0)

# static conditions at the fan face for the assumed Mach number
T_face = T0/(1+(gam-1)/2*M_face**2)
P_face = P0*(T_face/T0)**(gam/(gam-1))
rho_face = P_face/(R*T_face)
a_face = np.sqrt(gam*R*T_face)
V_face = M_face*a_face

m0 = rho_face*A*V_face

if __name__ == "__main__":
    print(f"Fan annulus area A = {A:.4f} m^2  (D_tip=1.549 m, htr=0.30)")
    print(f"Fan-face static: T={T_face:.1f} K, P={P_face/1000:.2f} kPa, rho={rho_face:.4f} kg/m^3, V={V_face:.1f} m/s")
    print(f"Preliminary mass flow estimate m0 = {m0:.1f} kg/s")