
import numpy as np
 
#  constants
R_air = 287.0
cp_c, gam_c = 1005.0, 1.40    		  # cold side (air)
cp_t, gam_t = 1244.0, 1.333    		 # hot side (combustion gas)
LHV = 43.1e6               	     		 # J/kg Jet-A
 
def atm(alt_m):
    T = 288.15 - 0.0065*alt_m
    P = 101325.0*(T/288.15)**5.2561
    return T, P
 
def Tt_Pt(T, P, M, gam):
    Tt = T*(1+(gam-1)/2*M**2)
    Pt = P*(1+(gam-1)/2*M**2)**(gam/(gam-1))
    return Tt, Pt
 
def crit_PR(gam):
    return ((gam+1)/2)**(gam/(gam-1))
 
def nozzle(Pt, Tt, P0, gam, R):
    PR, PRc = Pt/P0, crit_PR(gam)
    if PR > PRc:
        M, choked = 1.0, True
        T = Tt/(1+(gam-1)/2)
        P = Pt/PRc
    else:
        M, choked = np.sqrt(2/(gam-1)*(PR**((gam-1)/gam)-1)), False
        T, P = Tt/(1+(gam-1)/2*M**2), P0
    V = np.sqrt(gam*R*T)*M
    Veff = V + (R*T/V)*(1-P0/P) if choked else V
    return dict(M=M, T=T, P=P, V=V, Veff=Veff, choked=choked, PR=PR, PRc=PRc)
 
def cycle(flight, OPR, BPR, FPR, BoosterPR, T4, cool_frac=0.0,
          eta_f=0.88, eta_b_lpc=0.86, eta_c=0.87, eta_tH=0.90, eta_tL=0.91,
          eta_burn=0.99, eta_m=0.99, pi_d=0.995, pi_b=0.95, pi_n=0.98):
    """
    Real-cycle model with the booster (LPC) modelled as a separate stage on
    the LP spool, and a fraction of HPC exit air removed as turbine cooling
air """

    M0, alt = flight
    T0, P0 = atm(alt)
    a0 = np.sqrt(gam_c*R_air*T0)
    V0 = M0*a0
    Tt0, Pt0 = Tt_Pt(T0, P0, M0, gam_c)
 
    # inlet
    Pt2, Tt2 = Pt0*pi_d, Tt0
 
    # fan (bypass stream + core stream both pass through)
    Pt13 = Pt2*FPR
    Tt13s = Tt2*FPR**((gam_c-1)/gam_c)
    Tt13 = Tt2 + (Tt13s-Tt2)/eta_f
 
    # booster / LPC - core stream only, downstream of fan root
    Pt25 = Pt13*BoosterPR
    Tt25s = Tt13*BoosterPR**((gam_c-1)/gam_c)
    Tt25 = Tt13 + (Tt25s-Tt13)/eta_b_lpc

    # HPC
    HPC_PR = OPR/(FPR*BoosterPR)
    Pt3 = Pt25*HPC_PR
    Tt3s = Tt25*HPC_PR**((gam_c-1)/gam_c)
    Tt3 = Tt25 + (Tt3s-Tt25)/eta_c
    # burner (cooling air bled off here, bypasses combustion)
    Pt4 = Pt3*pi_b
    Tt4 = T4
    m_core_frac_burned = 1 - cool_frac       # fraction of core flow that is burned
    f = cp_t*(Tt4-Tt3) / (eta_burn*LHV - cp_t*Tt4)   # f defined per unit burned flow
 
    # HPT: drives HPC. Only the burned+fuel flow does turbine work; cooling
    # air is reintroduced at HPT exit (simple first-order treatment).
    w_hpc = cp_c*(Tt3-Tt25)                  # per unit TOTAL core flow
    dTt_hpt = (w_hpc/eta_m) / (m_core_frac_burned*(1+f)*cp_t)
    Tt45_gas = Tt4 - dTt_hpt                 # temperature of the burned gas leaving HPT rotor
    Tt45s = Tt4 - dTt_hpt/eta_tH
    Pt45 = Pt4*(Tt45s/Tt4)**(gam_t/(gam_t-1))
    # mix in cooling air (assumed to re-enter at Tt3, the HPC exit temperature)
    Tt45 = (m_core_frac_burned*(1+f)*cp_t*Tt45_gas + cool_frac*cp_c*Tt3) / \
           (m_core_frac_burned*(1+f)*cp_t + cool_frac*cp_c)
 
    # LPT: drives fan + booster
    w_fan = (1+BPR)*cp_c*(Tt13-Tt2)
    w_boost = cp_c*(Tt25-Tt13)
    dTt_lpt = ((w_fan+w_boost)/eta_m) / ((1+f)*cp_t)   # approx, full core mass now includes cooling air
    Tt5 = Tt45 - dTt_lpt
    Tt5s = Tt45 - dTt_lpt/eta_tL
    Pt5 = Pt45*(Tt5s/Tt45)**(gam_t/(gam_t-1))
 
    Pt9, Pt19 = Pt5*pi_n, Pt13*pi_n
    core = nozzle(Pt9, Tt5, P0, gam_t, R_air)
    byp  = nozzle(Pt19, Tt13, P0, gam_c, R_air)
 
    f_total = f*m_core_frac_burned             # fuel-air ratio referenced to total core flow
    Fsp = ((1+f_total)*core['Veff'] - V0 + BPR*(byp['Veff']-V0)) / (1+BPR)
    TSFC_SI = f_total/((1+BPR)*Fsp)
    TSFC_lb = TSFC_SI*2.20462*3600/0.224809
 
    return dict(T0=T0,P0=P0,V0=V0,Tt2=Tt2,Pt2=Pt2,Tt13=Tt13,Pt13=Pt13,
                Tt25=Tt25,Pt25=Pt25,Tt3=Tt3,Pt3=Pt3,Tt4=Tt4,Pt4=Pt4,f=f_total,
                Tt45=Tt45,Pt45=Pt45,Tt5=Tt5,Pt5=Pt5,Pt9=Pt9,Pt19=Pt19,
                core=core,byp=byp,Fsp=Fsp,TSFC_SI=TSFC_SI,TSFC_lb=TSFC_lb,HPC_PR=HPC_PR)
 
#  CFM56-7B26 baseline 
BASE = dict(OPR=32.6, BPR=5.1, FPR=1.60, BoosterPR=1.70, cool_frac=0.12)
 

cruise = cycle((0.80, 10668.0), T4=1480.0, **BASE)
takeoff = cycle((0.0, 0.0), T4=1630.0, **BASE)

print(" CRUISE M0.80 / 35,000ft ")
print("-----------------------------")
print(f"Fsp={cruise['Fsp']:.2f} N.s/kg  TSFC={cruise['TSFC_lb']:.4f} lb/(lbf.hr)  f={cruise['f']:.4f}")

print(" TAKE-OFF SL/M0 ")
print("-----------------------------")
print(f"Fsp={takeoff['Fsp']:.2f} N.s/kg  TSFC={takeoff['TSFC_lb']:.4f} lb/(lbf.hr)  f={takeoff['f']:.4f}")
 
F_target = 116988.0   # 26,300 lbf
m0_req = F_target/takeoff['Fsp']