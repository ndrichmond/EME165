import numpy as np
from matplotlib import pyplot as plt

#variable setup
T_air = 20 + 273.15 #kelvin
L_1 = 1.20 #meters, core length
L_2 = 0.25 #meters, handle length
r_core = 0.02 #meters, radius of core (and handle)
k_core = 0.35 #W/mK
k_handle = 0.05 #W/mK
R_tc = 10**-5 #m^2K/W, thermal resistance of handle-core boundary
R_tc_luke = 10**-4 #m^2K/W, thermal resistance of luke-handle boundary
T_luke = 37.1 + 273.15 #kelvin, temperature of Luke's hand
qdot_r0 = 10**6
h = 50 #W/m^2k
pi = np.pi

alpha = 0.8

def r(i): #radius at node i
    return (r_core/(nInd))*i

def qgen(r1,r2,dx): #calcluate the heat generated in a cylindrical volume
    return qdot_r0*pi*(dx)*((r2**2-r1**2) + ((1/(2*r_core**2))*(r1**4-r2**4)))

def aAnn(r1,r2): #calculate the area of an annulus
    return pi*(r2**2 - r1**2)

def aCyl(r,h):
    return 2*pi*r*h

#node variable setup
n = 5  #number of nodes in the radial direction
m = 200 #number of nodes in the axial direction in the core
M = 5  #number of nodes in the axial direction in the handle

nInd = n - 1 #for the purpose of indexing 
mInd = m - 1 #for the purpose of indexing 
MInd = M - 1 #for the purpose of indexing 

dr = r_core/(n-1)
dx_c = L_1/(m-1)
dx_h = L_2/(M-1)

T_core = np.zeros((n,m))
T_core_last = np.zeros((n,m))
t_core_updated = np.zeros((n,m))
T_handle = np.zeros((n,M))
T_handle_last = np.zeros((n,M))
t_handle_updated = np.zeros((n,M))

deltaT = 1 #maximum variation between iterations
tolerance = 0.01 #desired final tempreature variation between iterations

#iterations
iterationCounter = 0
while deltaT > tolerance and iterationCounter < 1000:
    #center edge nodes "corners"
    T_core[0,mInd] = (T_core_last[0,mInd-1]*(k_core*pi*(dr/2)**2/(dx_c/2)) + 
                      T_core_last[1,mInd]*(2*pi*k_core*dx_c/2) + 
                      qgen(0,dr/2,dx_c/2) + T_air*h*pi*(dr/2)**2*(dx_c/2)) / ((pi*(dr/2)**2*((2*k_core/2) + h)) + (2*pi*k_core*(dx_c/2)))
    T_core[0,0] = (T_core_last[0,1]*(k_core*aAnn(0,dr/2)/(dx_c/2)) + 
                   T_core_last[1,0]*(k_core*aCyl(dr/2,dx_c)/(dr/2)) + 
                   T_handle_last[0,0]*(aAnn(0,dr/2)/(R_tc)) + 
                   qgen(0,dr/2,dx_c/2)) / (pi*(dr/2)**2*((2*k_core/dx_c) + 1/R_tc) + pi*k_core*dx_c)
    
    #outer edge nodes "corners"
    T_core[nInd,mInd] = (T_core_last[nInd-1,mInd]*(k_core*aCyl((r_core - dr/2),dx_c/2)/(dr/2)) + 
                         T_core_last[nInd,mInd-1]*(k_core*aAnn((r_core - (dr/2)),r_core)/(dx_c/2)) + 
                         T_air*(h*(aCyl(r_core,(dx_c/2)) + aAnn(r_core,(r_core - (dr/2))))) + 
                         qgen((r_core - dr/2),r_core,dx_c/2)) / ((k_core*aCyl(r_core - (dr/2),dx_c/2))/(dr/2) + h*aCyl(r_core,dx_c/2) + pi*(r_core**2 - (r_core - (dr/2))**2)*(((2*k_core)/(dx_c)) + h))
    T_core[nInd,0] = (T_core_last[nInd-1,0]*(k_core*aCyl(r_core - (dr/2),dx_c/2)/(dr/2)) +
                      T_core_last[nInd,1]*(k_core*aAnn(r_core-(dr/2),r_core)/(dx_c/2)) +
                      T_handle_last[nInd,0]*(aAnn(r_core-(dr/2),r_core)/R_tc) +
                      T_air*h*aCyl(r_core,dx_c/2) + qgen(r_core - (dr/2),r_core,(dx_c/2))) / (k_core*aCyl(r_core - (dr/2),(dx_c/2)) + h*2*pi*r_core*(dx_c/2) + pi*(r_core**2 - (r_core - (dr/2))**2)*((k_core*dx_c/2) + (1/R_tc)))
    
    #core center nodes ("side")
    j = 1
    while j < mInd:
        T_core[0,j] = (T_core_last[1,j]*k_core*aCyl(dr/2,dx_c)/(dr/2) +
                       T_core_last[0,j-1]*k_core*aAnn(0,dr/2)/dx_c + 
                       T_core_last[0,j+1]*k_core*aAnn(0,dr/2)/dx_c + 
                       qgen(0,dr/2,dx_c)) / (k_core*aCyl(dr/2,dx_c)/(dr/2) + 2*k_core*aAnn(0,dr/2)/dx_c)
        j += 1

    #center of core handle contact nodes ("side")
    i = 1
    while i < nInd:
        T_core[i,0] = (T_core_last[i-1,0]*(k_core*aCyl(r(i) - (dr/2),dx_c/2)/dr) + 
                       T_core_last[i+1,0]*(k_core*aCyl(r(i) + (dr/2),dx_c/2)/dr) +
                       T_core_last[i,1]*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/(dx_c/2)) + 
                       T_handle_last[i,0]*(aAnn(r(i) - (dr/2),r(i) + (dr/2))/R_tc) + 
                       qgen(r(i) - (dr/2),r(i) + (dr/2),dx_c/2)) / ((k_core*aCyl(r(i) - (dr/2),dx_c/2)/dr) + (k_core*aCyl(r(i) + (dr/2),dx_c/2)/dr) + (k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/(dx_c/2)) + (aAnn(r(i) - (dr/2),r(i) + (dr/2))/R_tc))
        i += 1

    #edge of core nodes ("side")
    j = 1
    while j < mInd:
        T_core[nInd,j] = (T_core_last[nInd,j-1]*(k_core*aAnn(r(nInd) - (dr/2),r(nInd))/dx_c) + 
                          T_core_last[nInd,j+1]*(k_core*aAnn(r(nInd) - (dr/2),r(nInd))/dx_c) + 
                          T_core_last[nInd-1,j]*(k_core*aCyl(r(nInd) - (dr/2),dx_c)/(dr/2)) + 
                          T_air*(h*aCyl(r(nInd),dx_c)) +
                          qgen(r(nInd) - (dr/2),r(nInd),dx_c)) / ((k_core*aAnn(r(nInd) - (dr/2),r(nInd))/dx_c) + (k_core*aAnn(r(nInd) - (dr/2),r(nInd))/dx_c) + (k_core*aCyl(r(nInd) - (dr/2),dx_c)/(dr/2)) + (h*aCyl(r(nInd),dx_c)))
        j += 1

    #core tip nodes ("side")
    i = 1
    while i < nInd:
        T_core[i,mInd] = (T_core_last[i,mInd-1]*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/(dx_c/2)) +
                          T_core_last[i-1,mInd]*(k_core*aCyl(r(i) - (dr/2),dx_c/2)/dr) + 
                          T_core_last[i+1,mInd]*(k_core*aCyl(r(i) + (dr/2),dx_c/2)/dr) + 
                          T_air*(h*aAnn(r(i) - (dr/2),r(i) + (dr/2))) + 
                          qgen(r(i) - (dr/2),r(i) + (dr/2),dx_c/2)) / ((k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/(dx_c/2)) + (k_core*aCyl(r(i) - (dr/2),dx_c/2)/dr) + (k_core*aCyl(r(i) + (dr/2),dx_c/2)/dr) + (h*aAnn(r(i) - (dr/2),r(i) + (dr/2))))
        i += 1

    #interior nodes
    i = 1
    while i < nInd:
        j = 1
        while j < mInd:
            T_core[i,j] = (T_core_last[i-1,j]*(k_core*aCyl(r(i) - (dr/2),dx_c)/dr) + 
                           T_core_last[i+1,j]*(k_core*aCyl(r(i) + (dr/2),dx_c)/dr) + 
                           T_core_last[i,j-1]*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/dx_c) + 
                           T_core_last[i,j+1]*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/dx_c) + 
                           qgen(r(i) - (dr/2),r(i) + (dr/2),dx_c)) / ((k_core*aCyl(r(i) - (dr/2),dx_c)/dr) + (k_core*aCyl(r(i) + (dr/2),dx_c)/dr) + 2*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/dx_c))
            j += 1
        i += 1

    #the below is bugged somehow
    #repeat for the handle, but with these changes:
    # - instead of convection there is contact resistance between the handle and luke's skin
    # - use handle instead of core
    # - use MInd instead of mInd
    # - use dx_h instead of dx_c
    # - there is certainly a better way to do this with a function which i may implement

    #center edge nodes "corners"
    T_handle[0,MInd] = (T_handle_last[0,MInd-1]*(k_handle*aAnn(0,dr/2)/(dx_h/2)) + 
                        T_handle_last[1,MInd]*(k_handle*aCyl(dr/2,dx_h/2)/(dr/2)) + 
                        T_luke*(aAnn(0,dr/2)/R_tc_luke) +
                        qgen(0,dr/2,dx_h/2)) / ((k_handle*aAnn(0,dr/2)/(dx_h/2)) + (k_handle*aCyl(dr/2,dx_h/2)/(dr/2)) + (aAnn(0,dr/2)/R_tc_luke))

    T_handle[0,0] = (T_handle_last[0,1]*(k_handle*aAnn(0,dr/2)/(dx_h/2)) +
                     T_handle_last[1,0]*(k_handle*aCyl(dr/2,dx_h/2)/(dr/2)) + 
                     T_core_last[0,0]*(aAnn(0,dr/2)/R_tc) + 
                     qgen(0,dr/2,dx_h/2)) / ((k_handle*aAnn(0,dr/2)/(dx_h/2)) + (k_handle*aCyl(dr/2,dx_h/2)/(dr/2)) + (aAnn(0,dr/2)/R_tc))

    #outer edge nodes "coner"
    T_handle[nInd,MInd] = (T_handle_last[nInd-1,MInd]*(k_handle*aCyl(r_core - dr/2,dx_h/2)/(dr/2)) + 
                           T_handle_last[nInd,MInd-1]*(k_handle*aAnn(r_core - dr/2,r_core)/(dx_h/2)) + 
                           T_luke*((aCyl(r_core,dx_h/2) + aAnn(r_core + dr/2,r_core))/(R_tc_luke)) + 
                           qgen(r_core + dr/2,r_core,dx_h/2)) / ((k_handle*aCyl(r_core - dr/2,dx_h/2)/(dr/2)) + (k_handle*aAnn(r_core - dr/2,r_core)/(dx_h/2)) + ((aCyl(r_core,dx_h/2) + aAnn(r_core + dr/2,r_core))/(R_tc_luke)))
    
    T_handle[nInd,0] = (T_handle_last[nInd-1,0]*(k_handle*aCyl(r_core - dr/2,dx_h/2)/(dr/2)) + 
                        T_handle_last[nInd,1]*(k_handle*aAnn(r_core - dr/2,r_core)/(dx_h/2)) + 
                        T_luke*(aCyl(r_core,dx_h/2)/R_tc_luke) + 
                        T_core_last[nInd,0]*(aAnn(r_core - dr/2,r_core)/R_tc) + 
                        qgen(r_core - dr/2,r_core,dx_h/2)) / ((k_handle*aCyl(r_core - dr/2,dx_h/2)/(dr/2)) + (k_handle*aAnn(r_core - dr/2,r_core)/(dx_h/2)) + (aCyl(r_core,dx_h/2)/R_tc_luke) + (aAnn(r_core - dr/2,r_core)/R_tc))
    
    #t_core_updated = T_core.copy()
    #t_handle_updated = T_handle.copy()

    deltaT = np.amax(np.absolute(T_core_last - T_core))

    T_core_last = T_core.copy()
    #T_handle_last = alpha*t_handle_updated.copy() + (1-alpha) * T_handle_last.copy()

    iterationCounter += 1
    print(iterationCounter,deltaT)
    #print(T_core)


""" 
    #center edge nodes "corners"
    T_core[0,MInd] = (T_core_last[0,mInd-1]*(k_core*pi*(dr/2)**2/(dx_c/2)) + 
                      T_core_last[1,mInd]*(2*pi*k_core*dx_c/2) + 
                      qgen(0,dr/2,dx_c/2) + h*pi*(dr/2)**2*(dx_c/2)) / ((pi*(dr/2)**2*((2*k_core/2) + h)) + (2*pi*k_core*(dx_c/2)))
    T_core[0,0] = (T_core_last[0,1]*(k_core*aAnn(0,dr/2)/(dx_c/2)) + 
                   T_core_last[1,0]*(k_core*aCyl(dr/2,dx_c)/(dr/2)) + 
                   T_handle_last[0,0]*(aAnn(0,dr/2)/(R_tc)) + 
                   qgen(0,dr/2,dx_c/2)) / (pi*(dr/2)**2*((2*k_core/dx_c) + 1/R_tc) + pi*k_core*dx_c)
    
    #outer edge nodes "corners"
    T_core[nInd,mInd] = (T_core_last[nInd-1,mInd]*(k_core*aCyl((r_core - dr/2),dx_c/2)/(dr/2)) + 
                         T_core_last[nInd,mInd-1]*(k_core*aAnn((r_core - (dr/2)),r_core)/(dx_c/2)) + 
                         T_air*(h*(aCyl(r_core,(dx_c/2)) + aAnn(r_core,(r_core - (dr/2))))) + 
                         qgen((r_core - dr/2),r_core,dx_c/2)) / ((k_core*aCyl(r_core - (dr/2),dx_c/2))/(dr/2) + h*aCyl(r_core,dx_c/2) + pi*(r_core**2 - (r_core - (dr/2))**2)*(((2*k_core)/(dx_c)) + h))
    T_core[nInd,0] = (T_core_last[nInd-1,0]*(k_core*aCyl(r_core - (dr/2),dx_c/2)/(dr/2)) +
                      T_core_last[nInd,1]*(k_core*aAnn(r_core-(dr/2),r_core)/(dx_c/2)) +
                      T_handle_last[nInd,0]*(aAnn(r_core-(dr/2),r_core)/R_tc) +
                      T_air*h*aCyl(r_core,dx_c/2) + qgen(r_core - (dr/2),r_core,(dx_c/2))) / (k_core*aCyl(r_core - (dr/2),(dx_c/2)) + h*2*pi*r_core*(dx_c/2) + pi*(r_core**2 - (r_core - (dr/2))**2)*((k_core*dx_c/2) + (1/R_tc)))
    
    #core center nodes ("side")
    j = 1
    while j < mInd:
        T_core[0,j] = (T_core_last[1,j]*k_core*aCyl(dr/2,dx_c)/(dr/2) +
                       T_core_last[0,j-1]*k_core*aAnn(0,dr/2)/dx_c + 
                       T_core_last[0,j+1]*k_core*aAnn(0,dr/2)/dx_c + 
                       qgen(0,dr/2,dx_c)) / (k_core*aCyl(dr/2,dx_c)/(dr/2) + 2*k_core*aAnn(0,dr/2)/dx_c)
        j += 1

    #center of core handle contact nodes ("side")
    i = 1
    while i < nInd:
        T_core[i,0] = (T_core_last[i-1,0]*(k_core*aCyl(r(i) - (dr/2),dx_c/2)/dr) + 
                       T_core_last[i+1,0]*(k_core*aCyl(r(i) + (dr/2),dx_c/2)/dr) +
                       T_core_last[i,1]*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/(dx_c/2)) + 
                       T_handle_last[i,0]*(aAnn(r(i) - (dr/2),r(i) + (dr/2))/R_tc) + 
                       qgen(r(i) - (dr/2),r(i) + (dr/2),dx_c/2)) / ((k_core*aCyl(r(i) - (dr/2),dx_c/2)/dr) + (k_core*aCyl(r(i) + (dr/2),dx_c/2)/dr) + (k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/(dx_c/2)) + (aAnn(r(i) - (dr/2),r(i) + (dr/2))/R_tc))
        i += 1

    #edge of core nodes ("side")
    j = 1
    while j < mInd:
        T_core[nInd,j] = (T_core_last[nInd,j-1]*(k_core*aAnn(r(nInd) - (dr/2),r(nInd))/dx_c) + 
                          T_core_last[nInd,j+1]*(k_core*aAnn(r(nInd) - (dr/2),r(nInd))/dx_c) + 
                          T_core_last[nInd-1,j]*(k_core*aCyl(r(nInd) - (dr/2),dx_c)/(dr/2)) + 
                          T_air*(h*aCyl(r(nInd),dx_c)) +
                          qgen(r(nInd) - (dr/2),r(nInd),dx_c)) / ((k_core*aAnn(r(nInd) - (dr/2),r(nInd))/dx_c) + (k_core*aAnn(r(nInd) - (dr/2),r(nInd))/dx_c) + (k_core*aCyl(r(nInd) - (dr/2),dx_c)/(dr/2)) + (h*aCyl(r(nInd),dx_c)))
        j += 1

    #core tip nodes ("side")
    i = 1
    while i < nInd:
        T_core[i,mInd] = (T_core_last[i,mInd-1]*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/(dx_c/2)) +
                          T_core_last[i-1,mInd]*(k_core*aCyl(r(i) - (dr/2),dx_c/2)/dr) + 
                          T_core_last[i+1,mInd]*(k_core*aCyl(r(i) + (dr/2),dx_c/2)/dr) + 
                          T_air*(h*aAnn(r(i) - (dr/2),r(i) + (dr/2))) + 
                          qgen(r(i) - (dr/2),r(i) + (dr/2),dx_c/2)) / ((k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/(dx_c/2)) + (k_core*aCyl(r(i) - (dr/2),dx_c/2)/dr) + (k_core*aCyl(r(i) + (dr/2),dx_c/2)/dr) + (h*aAnn(r(i) - (dr/2),r(i) + (dr/2))))
        i += 1

    #interior nodes
    i = 1
    while i < nInd:
        j = 1
        while j < mInd:
            T_core[i,j] = (T_core_last[i-1,j]*(k_core*aCyl(r(i) - (dr/2),dx_c)/dr) + 
                           T_core_last[i+1,j]*(k_core*aCyl(r(i) + (dr/2),dx_c)/dr) + 
                           T_core_last[i,j-1]*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/dx_c) + 
                           T_core_last[i,j+1]*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/dx_c) + 
                           qgen(r(i) - (dr/2),r(i) + (dr/2),dx_c)) / ((k_core*aCyl(r(i) - (dr/2),dx_c)/dr) + (k_core*aCyl(r(i) + (dr/2),dx_c)/dr) + 2*(k_core*aAnn(r(i) - (dr/2),r(i) + (dr/2))/dx_c))
            j += 1
        i += 1
"""

#plot everything

T_core = T_core.copy() - 273

Nr, Nz = T_core.shape

r = np.linspace(0, r_core, Nr)
z = np.linspace(0, L_1, Nz)
R, Z = np.meshgrid(r, z, indexing='ij')  # use 'ij' for (r, z) ordering

T_mirrored = np.vstack([np.flipud(T_core), T_core])  # shape: (2*Nr, Nz)
r_mirrored = np.hstack([-np.flip(r), r])   # shape: (2*Nr,)
R_mirrored, Z_mirrored = np.meshgrid(r_mirrored, z, indexing='ij')

# --- Plot ---
plt.figure(figsize=(8, 5))
pcm = plt.pcolormesh(Z_mirrored, R_mirrored, T_mirrored, cmap='coolwarm', shading='auto')
plt.colorbar(pcm, label='Temperature (°C)')
plt.xlabel('Axial position z (m)')
plt.ylabel('Radial position r (m)')
plt.title('Mirrored Temperature Distribution in Cylinder')
plt.axis('equal')
plt.tight_layout()
plt.show()