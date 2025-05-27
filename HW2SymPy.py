import sympy as sym
import numpy as np

Jp = sym.Symbol("Jp")
Jf = sym.Symbol("Jf")
Jw = sym.Symbol("Jw")

e_p = 0.8
e_f = 0.9

T_p = 300
T_f = 800

Ri = 0.15/0.3
Rj = Ri

S = 1 + (1+Rj**2)/Ri**2

Fij = (1/2)*(S-(S**2-4)**(1/2))

Fpw = 1-Fij
Ffw = 1-Fij
Fpf = Fij
Ffp = Fij



sigma = 5.67*10**(-8)

Ebp = sigma*T_p**4
Ebf = sigma*T_f**4

Ap = np.pi*(0.15**2)
Af = Ap

solutions = sym.solve([(Ebp - Jp)/((1-e_p)/(e_p*Ap))  - (((Jp - Jw)/((1/(Ap*Fpw)))) + ((Jp-Jf)/(1/(Ap*Fpf)))),
                 (Ebf - Jf)/((1-e_f)/(e_f*Af)) - (((Jf - Jp)/(1/(Af*Ffp))) + ((Jf-Jw)/(1/(Af*Ffw)))),
                 (Ebp - Jp)/((1-e_p)/(e_p*Ap))  + ((Ebf-Jf)/((1-e_f)/(e_f*Af)))],[Jp,Jf,Jw])

Jp = solutions[Jp]
Jw = solutions[Jw]
#print(Jw)
#print((Ebp - Jp)/((1-e_p)/(e_p*Ap)))

Tw = (Jw/sigma)**(1/4)
#print(Tw)

#print(2*np.pi*0.15*0.3)
#print(np.pi*0.15**2)

Gf = sym.Symbol("Gf")
Gp = sym.Symbol("Gp")

Fwp = Fpw / (2*np.pi*0.15*0.3/(np.pi*0.15**2))
Fwf = Fwp

solns = sym.solve    ( [Fwp*(sigma*Tw**4) + Ffp*(e_f*sigma*T_f**4 + (1 - e_f)*Gf) - Gp,
                     Fwf*(sigma*Tw**4) + Fpf*(e_p*sigma*T_p**4 + (1 - e_p)*Gp) - Gf], 
            [Gp,Gf])
Gf = solns[Gf]
Gp = solns[Gp]

Gw = sigma*Tw**4

print(Gf*Ffp + Fwp*Gw)

G = (Jp - e_p*sigma*T_p**4)/(1-e_p)
print(G)
