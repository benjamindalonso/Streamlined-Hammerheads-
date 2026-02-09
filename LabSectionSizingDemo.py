import numpy as np
import matplotlib.pyplot as plt

# Parameter definitions (Add any additional variables to this list)

AR = 2.5        # From openVSP
s = 40          # Wing span from openVSP
s_ref = 740     # From openVSP
S_wet = 3144    # Estimated wetted area from DragPolars.py
c_f = 0.0040    # From table on resources doc
C_D_0 = c_f * (S_wet / s_ref)
print("Zero-lift drag coefficient C_D_0:", C_D_0)
N_eng = 1  # Number of engines
k_s = 1.1  # Based on RFP (10% above stall speed)
C_L_max = 1.6 # An estimate based on other fighter aircraft
G = 0.15  # Gradient (%) for fighter jet 
e = 0.75  # Oswald efficiency factor for fighter jet
WS = np.linspace(1,100,100)
W0 = 67822  # Maximum takeoff weight from Assignment 2
rhoRatio = 1 # assume standard day at sea level
TOP25 = 1092/37.5 # Parametric value for takeoff performance on a cvn-78 carrier 
Sland = 4000 # Target landing distance in feet without arresetment 
g = 32.17
rho_turn = 0.001267  # air density at 20,000 ft (slugs/ft^3)
V_turn = 500        # Assume 500 ft/s from PrelimSizing code
q = 0.5 * rho_turn * (V_turn**2) # Dynamic pressure
Turn_rate = 0.1745 # in Rad/s from rfp's 10 deg/s
k = 1 / (np.pi * e * AR) # Induced drag factor
n = np.sqrt(1 + ((V_turn * Turn_rate) / g)**2)



# These equations need to be filled in after you find out the relationship between T/w and W/S for each constraint
TW_takeoff = WS / (TOP25 * rhoRatio * C_L_max) # done by Mark
WS_landing = ((rhoRatio * C_L_max)/(80*0.6)) * (Sland-1000)
# Climb Constraint ======================================================
k = 1/(np.pi*AR*e)
def climb_constraint(MTOW,C_D0,CLmax,MLW,airdens,k,vert_speed):
    k_s = 1.1  # stall safety factor
    # C_D0 is zero drag lift
    C_L_TO = CLmax # take-off maximum lift
    W = MTOW
    W_land = MLW # mean landing weight
    # airden is air density
    # k is not defined in slides
    # vert_speed is in ft/s
    G = (vert_speed)/(((2*k_s*W)/(airdens*C_L_TO)-(vert_speed)**2)**(1/2)) # Climb Gradient
    Climb_TW = ((k_s**2)*C_D0)/C_L_TO + (k*C_L_TO)/(k_s**2) + G # Climb thrust to weight
    Climb_TW = Climb_TW*(W_land/W)*(1/0.80) # climb thrust to weight with correction factors
    return Climb_TW
#========================================================================
TW_climb = WS*climb_constraint(67822,C_D_0,C_L_max,W0,0.02377,k,5000/60)/WS
TW_turn = q * (C_D_0 / WS + k * (n / q)**2 * WS) #done by Aiden

# ===================== cruise constraint (M=1.6 30,000 ft) =====================

def cruise(alt_ft):
    g0 = 32.174       # ft/s^2
    R  = 1716.59      # ft*lbf/(slug*R)
    T0 = 518.67       # R
    p0 = 2116.22      # lbf/ft^2
    L  = -0.00356616  # R/ft

    T = T0 + L * alt_ft
    p = p0 * (T / T0) ** (-g0 / (L * R))
    rho = p / (R * T)
    return T, rho

M_cruise = 1.6
alt_cruise_ft = 30000.0

T_cr, rho_cr = cruise(alt_cruise_ft)
gamma = 1.4
a_cr = np.sqrt(gamma * 1716.59 * T_cr)  # ft/s
V_cr = M_cruise * a_cr                  # ft/s
q_cr = 0.5 * rho_cr * V_cr**2           # lbf/ft^2
k_ind = 1.0 / (np.pi * e * AR)
beta_cruise  = 0.90  
alpha_cruise = 0.75   
WS_cr = beta_cruise * WS
TW_cr = (q_cr * C_D_0) / WS_cr + (k_ind * WS_cr) / q_cr
TW_cruise = (beta_cruise / alpha_cruise) * TW_cr

# ================================================================================

plt.figure(figsize=(8,4))
plt.title('T/W - W/S')
plt.xlabel("W/S $(lb/ft^2)$")
plt.ylabel("T/W")
plt.plot(WS, TW_takeoff, label='Takeoff field length', linestyle='-', linewidth=2)
plt.axvline(x=WS_landing, color='green', linestyle='--', linewidth=2, label=f'Landing')
plt.plot(WS,TW_climb,color = 'blue', linestyle='-.',linewidth=2, label='200 FPM Initial Climb')
plt.plot(WS, TW_cruise, color='red', linestyle=':', linewidth=2, label='Cruise (M=1.6 @ 30kft)')
plt.plot(WS,TW_turn,color = 'orange', linestyle='-.',linewidth=2, label='Turn Constraint')
plt.ylim(0, 1.5)
plt.legend(loc='best')
plt.show()