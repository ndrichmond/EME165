import sympy as sym
import numpy as np

sigma = 5.67e-8

alpha = 0.96
h_i = 20
h_o = 45
r_i = 4
t_w = 0.01
t_ice = 0.1
k_w = 0.5
k_ice = 2.22
R_tc = 5*10**(-2)

'''
alpha = sym.Symbol("alpha")
h_i = sym.Symbol("h_i")
h_o = sym.Symbol("h_o")
r_i = sym.Symbol("r_i")
t_w = sym.Symbol("t_w")
t_ice = sym.Symbol("t_ice")
k_w = sym.Symbol("k_w")
k_ice = sym.Symbol("k_ice")
R_tc = sym.Symbol("R_tc")
'''

T_sky = 258
T_info = 265
#T_sky = sym.Symbol("T_sky")
#T_info = sym.Symbol("T_info")
T_s_o = sym.Symbol("T_s_o",real=True,positive=True)


A_wi = (1/2)*4*np.pi*r_i**2
A_ww = (1/2)*4*np.pi*(r_i + t_w)**2
A_o = (1/2)*4*np.pi*(r_i + t_w + t_ice)**2

#A_wi = (1/2)*4*pi*r_i**2
#A_ww = (1/2)*4*pi*(r_i + t_w)**2
#A_o = (1/2)*4*pi*(r_i + t_w + t_ice)**2

R_total = (1/(h_i*A_wi)) + (((1 - ((r_i)/(r_i+t_w)))/(2*np.pi*k_w*r_i))) + R_tc/A_ww + (1-((r_i+t_w)/(r_i + t_w + t_ice)))/(2*np.pi*k_ice*(r_i+t_w))
#print(R_total)
print(R_total)
T_infi = 298

solns = sym.solve([alpha*sigma*T_s_o**4 + h_o*A_o*(T_s_o-T_info) - alpha*sigma*T_sky**4 - ((T_s_o - T_infi)/(R_total))])
T_s_o = solns[0][T_s_o]
print(T_s_o)
q_wall = alpha*sigma*T_s_o**4 + h_o*A_o*(T_s_o-T_info) - alpha*sigma*T_sky**4
print(q_wall)
#for floor q calculations

T_soil = 270
T_floor = sym.Symbol("T_floor")
t_f_ice = 5
A_f = np.pi*r_i**2

q_floor = (T_infi - T_soil)/((1/(h_i*A_f)) + (t_w/(k_w*A_f)) + (R_tc/A_f) + (t_f_ice/(k_ice*A_f)))
print(q_floor)


