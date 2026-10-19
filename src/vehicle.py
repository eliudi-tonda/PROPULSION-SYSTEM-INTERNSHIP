# variable-property cycle

"""Independent variable-property (T-dependent cp) check of the SLS take-off cycle.
cp_air(T) polynomial (Walsh & Fletcher form), combustion-gas correction proportional to FAR.
Same component efficiencies/losses as cycle.py so the ONLY change is the gas model."""

import numpy as np
from scipy.optimize import brentq

A=[0.992313,0.236688,-1.852148,6.083152,-8.893933,7.097112,-3.234725,0.794571,-0.081873]
B=[-0.718874,8.747481,-15.863157,17.254096,-10.233795,3.081778,-0.361112,-0.003919]
R=287.0; TREF=298.15; LHV=43.1e6

def _hp(c,th): return sum(c[i]*th**(i+1)/(i+1) for i in range(len(c)))

def _sp(c,th): return c[0]*np.log(th)+sum(c[i]*th**i/i for i in range(1,len(c)))

def cp(T,f=0.0):
    th=T/1000.0
    return 1000*(sum(A[i]*th**i for i in range(len(A)))+ f/(1+f)*sum(B[i]*th**i for i in range(len(B))))

def h(T,f=0.0):   # J/kg relative to TREF
    th,tr=T/1000.0,TREF/1000.0
    ha=1e6*(_hp(A,th)-_hp(A,tr)); hb=1e6*(_hp(B,th)-_hp(B,tr))
    return ha+f/(1+f)*hb

def phi(T,f=0.0):
    th=T/1000.0
    return 1000*(_sp(A,th)+f/(1+f)*_sp(B,th))

def T_from_h(hv,f=0.0): return brentq(lambda T:h(T,f)-hv,150,3000)

def comp(Tin,PR,eta,f=0.0):
    Ts=brentq(lambda T:phi(T,f)-phi(Tin,f)-R*np.log(PR),Tin,3000)
    hout=h(Tin,f)+(h(Ts,f)-h(Tin,f))/eta
    return T_from_h(hout,f)


def turb_drop(Tin,dh,eta,f,Pin):
    hs=h(Tin,f)-dh/eta
    Ts=T_from_h(hs,f); Pout=Pin*np.exp((phi(Ts,f)-phi(Tin,f))/R)
    Tout=T_from_h(h(Tin,f)-dh,f); return Tout,Pout

def nozzle(Pt,Tt,P0,f):
    # throat (M=1)
    g=lambda T: cp(T,f)/(cp(T,f)-R)
    Tth=brentq(lambda T:h(Tt,f)-h(T,f)-0.5*g(T)*R*T,200,Tt)
    Pth=Pt*np.exp((phi(Tth,f)-phi(Tt,f))/R)
    if Pth>P0:
        V=np.sqrt(g(Tth)*R*Tth); return V+(R*Tth/V)*(1-P0/Pth),True,Tth,Pth
    Te=brentq(lambda T:phi(T,f)-phi(Tt,f)-R*np.log(P0/Pt),200,Tt)
    return np.sqrt(2*(h(Tt,f)-h(Te,f))),False,Te,P0

def cycle(T4,OPR=27.2,BPR=5.1,FPR=1.60,BPRb=1.70,cool=0.12,eta_f=0.88,eta_b=0.86,eta_c=0.87,eta_tH=0.90,eta_tL=0.91,
          eta_burn=0.99,eta_m=0.99,pi_d=0.995,pi_b=0.95,pi_n=0.98,P0=101325.0,T0=288.15):
    Tt2=T0;Pt2=P0*pi_d
    Tt13=comp(Tt2,FPR,eta_f);Pt13=Pt2*FPR
    Tt25=comp(Tt13,BPRb,eta_b);Pt25=Pt13*BPRb
    HPC=OPR/(FPR*BPRb);Tt3=comp(Tt25,HPC,eta_c);Pt3=Pt25*HPC;Pt4=Pt3*pi_b
    # fuel-air (burner) : h_a(Tt3)+f*eta*LHV=(1+f)*h_g(T4,f)
    fb=brentq(lambda f:h(Tt3,0)+f*eta_burn*LHV-(1+f)*h(T4,f),0.001,0.06)
    mb=1-cool; ftot=fb*mb
    whpc=h(Tt3)-h(Tt25)
    dh_hpt=(whpc/eta_m)/(mb*(1+fb))
    Tt45g,Pt45=turb_drop(T4,dh_hpt,eta_tH,fb,Pt4)
    hmix=(mb*(1+fb)*h(Tt45g,fb)+cool*h(Tt3,0))/(mb*(1+fb)+cool)
    # composition of mixed stream
    Tt45=T_from_h(hmix,ftot)
    wlp=((1+BPR)*(h(Tt13)-h(Tt2))+(h(Tt25)-h(Tt13)))/eta_m
    dh_lpt=wlp/(1+ftot)
    Tt5,Pt5=turb_drop(Tt45,dh_lpt,eta_tL,ftot,Pt4*0+Pt45*1.0)
    Ve9,ch9,T9,P9=nozzle(Pt5*pi_n,Tt5,P0,ftot)
    Ve19,ch19,T19,P19=nozzle(Pt13*pi_n,Tt13,P0,0.0)
    Fsp=((1+ftot)*Ve9+BPR*Ve19)/(1+BPR)
    TSFC=ftot/((1+BPR)*Fsp)*2.20462*3600/0.224809
    return dict(Tt2=Tt2,Pt2=Pt2,Tt13=Tt13,Pt13=Pt13,Tt25=Tt25,Pt25=Pt25,Tt3=Tt3,Pt3=Pt3,Tt4=T4,Pt4=Pt4,Tt45=Tt45,Pt45=Pt45,
                Tt5=Tt5,Pt5=Pt5,Pt9=Pt5*pi_n,Pt19=Pt13*pi_n,T9=T9,P9=P9,T19=T19,Tt45g=Tt45g,
                f=ftot,fb=fb,Fsp=Fsp,TSFC=TSFC,ch9=ch9,ch19=ch19,Ve9=Ve9,Ve19=Ve19,cool=cool,BPR=BPR,HPC_PR=HPC)

if __name__=="__main__":
    import cycle as C
    for (T4,OPR) in ((1615,32.6),(1635,27.2)):
        v=cycle(T4,OPR=OPR); p=dict(C.BASE);p['OPR']=OPR;p.update(C.REF)
        c=C.cycle((0,0),T4=T4,**p)
        print(T4,OPR,'VAR Fsp %.2f TSFC %.4f f %.4f Tt3 %.1f Tt45 %.1f Tt5 %.1f'%(v['Fsp'],v['TSFC'],v['f'],v['Tt3'],v['Tt45'],v['Tt5']),v['ch9'],v['ch19'])
        print(T4,OPR,'CST Fsp %.2f TSFC %.4f f %.4f Tt3 %.1f Tt45 %.1f Tt5 %.1f'%(c['Fsp'],c['TSFC_lb'],c['f'],c['Tt3'],c['Tt45'],c['Tt5']))