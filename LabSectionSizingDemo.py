import numpy as np
import matplotlib.pyplot as plt

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
G = 0.012  # Gradient (%)
e = 0.8  # Oswald efficiency factor
coef_1_climb = (1/0.8) * (N_eng / (N_eng - 1)) * ((k_s**2) / C_L_max * C_D_0 + C_L_max / (np.pi * AR * e * k_s**2) + G)
print("Climb gradient coefficient:", coef_1_climb)

WS = np.linspace(1,300,30)

TW_takeoff = 0.001644*WS
TW_landing = 295.75*np.ones(30)
TW_climb = coef_1_climb*np.ones(30)
TW_cruise = 8.4566/WS + 0.0003865*WS


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

