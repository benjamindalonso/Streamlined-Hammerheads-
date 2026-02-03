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
coef_1_climb = (1/e) * (N_eng / (N_eng - 1)) * ((k_s**2) / C_L_max * C_D_0 + C_L_max / (np.pi * AR * e * k_s**2) + G)
print("Climb gradient coefficient:", coef_1_climb)
WS = np.linspace(1,300,30)
W0 = 67822  # Maximum takeoff weight from Assignment 2
rhoRatio = 1 # assume standard day at sea level
TOP25 = 1092/37.5 # Parametric value for takeoff performance on a cvn-78 carrier 


# These equations need to be filled in after you find out the relationship between T/w and W/S for each constraint
TW_takeoff = WS / (TOP25 * rhoRatio * C_L_max)
TW_landing = 
TW_climb = 
TW_cruise = 


plt.figure(figsize=(8,4))
plt.title('T/W - W/S')
plt.xlabel("W/S $(lb/ft^2)$")
plt.ylabel("T/W")
plt.plot(WS, TW_takeoff, label='Takeoff field length', linestyle='-', linewidth=2)
plt.plot(TW_landing, np.linspace(0,1,30), label='Landing field length', linestyle='-', linewidth=2)
plt.plot(WS, TW_climb, label='Takeoff climb', linestyle='-', linewidth=2)
plt.plot(WS, TW_cruise, label='Cruise', linestyle='-', linewidth=2)
plt.ylim(0, 0.5)
plt.legend(loc='best')
plt.show()

