from cmath import exp

# Crew and Payload Weights

Wcrew = 400          # Crew weight in lbs
Wpayload = 6802      # Payload weight in lbs (based on RFP requirement for armament)


# Regression constants

A = 2.392

c = -.13

# Mission segment weight fractions

warmup   = 0.99

taxi     = 0.99

takeoff  = 0.99

climb    = 0.96             

combat   = 0.96             # Assumed same as climb

descent  = 0.99

landing  = 0.995

goaround = 0.99             # This would be if landing is aborted and pilot has to re attempt


# Breguet Equations (Cruise and Loiter)

import numpy as np

Cl_Cd = 14   # Lift-to-Drag Ratio during cruise and loiter (14 for design B, 9 for A)
R = 2000            # Range in nmi 
ct = 0.8            # Thrust specific fuel consumtion in nmi/hr (0.8 for design B, 1 for A)
V = 550             # Cruise Velocity in nmi/hr (550 for design B, 1000 for A)
E = .333            # Time Spent Loitering in hours

cruise = np.exp((-R*ct) / (V*Cl_Cd))  

loiter = np.exp((-E*ct) / (Cl_Cd))    

# Mission weight fraction

Wn_W0 = (warmup * taxi * takeoff * climb * cruise * combat * loiter * descent * goaround * goaround * landing)       # 2 goaround allowance per the RFP requirements

# Fuel fractions

F_used = 1 - Wn_W0          # Used fuel fraction
F = 1.06 * F_used           # Total fuel fraction with reserves

# Initial Weight Guess and Convergence Loop

Wo = 80000           # Initial guess for Takeoff Gross Weight in lbs
Wo_history = []             # To store Wo values for convergence plot
previous_Wo = 0             # To track previous Wo for convergence check
err = 1e-6  
delta = 2*err                # Convergence error tolerance

while delta > err:              # While loop to solve for takeoff weight
    Wo_history.append(Wo)

    We_Wo = A * Wo ** c
    Wo_new = (Wcrew + Wpayload) / (1 - F - We_Wo)
    delta = abs(Wo_new - Wo) / abs(Wo_new)
    Wo = Wo_new

Wo_history = np.array(Wo_history)


# Final weights

We = We_Wo * Wo      # Empty Weight
Wfuel_total = F * Wo          # Total fuel weight
Wfuel_used = F_used * Wo      # Used Fuel weight
Wfuel_reserved = Wfuel_total - Wfuel_used      # Reserve Fuel weight
W_landing = Wn_W0 * Wo        # Landing weight

empty_weight_fraction = We / Wo        # Empty weight fraction
empty_weight_fraction_percent = empty_weight_fraction * 100  # Empty weight fraction in percent

# Weight Estimates Output
print("Lift to Drag Ratio (Cl/Cd): " + str(Cl_Cd))
print(f"\n--- Weight Estimates ---") # Print Weight Estimates Header
print("Gross Takeoff Weight: " + str(round(Wo)) + " lbs")
print(f"Empty Weight: {We:.2f} lbs") # Print Empty Weight
print("Regression's Empty Weight Fraction (We/W0): " + str(round(empty_weight_fraction ,3)))
print(f"Regression's Empty Weight Percentage: {empty_weight_fraction_percent:.2f}%") # Print Regression's Empty Weight Percentage
print(f"Landing Weight: {W_landing:.2f} lbs") # Print Landing Weight
print(f"Total Fuel Weight: {Wfuel_total:.2f} lbs") # Print Total Fuel Weight
print(f"Used Fuel Weight: {Wfuel_used:.2f} lbs") # Print Used Fuel Weight
print(f"Reserve Fuel Weight: {Wfuel_reserved:.2f} lbs") # Print Reserve Fuel Weight
print("Crew Weight: " + str(round(Wcrew)) + " lbs")
print("Payload Weight: " + str(round(Wpayload)) + " lbs")
print("Warm up Fuel Fraction (Wn/W0): " + str(round(warmup, 3)))
print("Taxi Fuel Fraction (Wn/W0): " + str(round(taxi, 3)))
print("Takeoff Fuel Fraction (Wn/W0): " + str(round(takeoff, 3)))
print("Cruise Fuel Fraction (Wn/W0): " + str(round(cruise, 3)))
print("Combat Fuel Fraction (Wn/W0): " + str(round(combat, 3)))
print("Loiter Fuel Fraction (Wn/W0): " + str(round(loiter, 3)))
print("Go Around Fuel Fraction (Wn/W0): " + str(round(goaround, 3)))
print("Landing Fuel Fraction (Wn/W0): " + str(round(landing, 3)))

# Plot Convergence
import matplotlib.pyplot as plt
# Plot convergence
plt.figure(figsize=(8, 4))
plt.title("Takeoff Weight Convergence")
plt.xlabel("Iteration")
plt.ylabel("W₀ (lbf)")
plt.plot(Wo_history, label="W₀", linewidth=2)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()



