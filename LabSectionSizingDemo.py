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


# These equations need to be filled in after you find out the relationship between T/w and W/S for each constraint
TW_takeoff = WS / (TOP25 * rhoRatio * C_L_max) # done by Mark
WS_landing = ((rhoRatio * C_L_max)/(80*0.6)) * (Sland-1000)
TW_climb = 0.1 # Change after you find this relationship
TW_cruise = 0.1 # Change after you find this relationship


plt.figure(figsize=(8,4))
plt.title('T/W - W/S')
plt.xlabel("W/S $(lb/ft^2)$")
plt.ylabel("T/W")
plt.plot(WS, TW_takeoff, label='Takeoff field length', linestyle='-', linewidth=2)
plt.axvline(x=WS_landing, color='green', linestyle='--', linewidth=2, label=f'Landing')
# Add your other contraints here
plt.ylim(0, 1.5)
plt.legend(loc='best')
plt.show()

