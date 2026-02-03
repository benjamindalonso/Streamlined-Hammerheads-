import numpy as np
import matplotlib.pyplot as plt

# Aircraft parameters
Wo = 73993 # Initial takeoff weight in lbs (From weight estimation code)
Sinitial = 740 # Initial wing area in ft²
WS_actual = Wo / Sinitial # Actual wing loading in lb/ft²

# Design assumptions
rho_ratio = 1.0 # Assume density of takeoff is sea level standard
TOP25 = (1092/37.5) # Parametric value for takeoff performance based on CVN 78 runway length
CLmax = 1.6  # Max lift coefficient during takeoff (Assumption)

# W/S range
WS = np.linspace(20, 200, 100)

# Actual takeoff line
TW_actual = WS / (rho_ratio * CLmax * TOP25)

plt.figure()
plt.plot(WS, TW_actual, 'b', label="Takeoff")

plt.xlabel("W/S (lb/ft²)")
plt.ylabel("T/W")
plt.title("Takeoff Constraint")
plt.legend()
plt.grid(True)
plt.show()

print("W/S (lb/ft²):", WS_actual)
print("T/W:", TW_actual[np.argmin(np.abs(WS - WS_actual))])