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

climb    = 0.94

combat   = 0.92

descent  = 0.99

landing  = 0.995

goaround = 0.99             # This would be if landing is aborted and pilot has to re attempt


# Breguet Equations (Cruise and Loiter)

import numpy as np

Cl_Cd = 14          # Lift-to-Drag Ratio during cruise and loiter
R = 2000            # Range in nmi
ct = 0.8            # Thrust specific fuel consumtion in nmi/hr
V = 550             # Velocity in nmi/hr
E = .333            # Time Spent Loitering in hours

cruise = np.exp((-R*ct) / (V*Cl_Cd))  

loiter = np.exp((-E*ct) / (Cl_Cd))    

# Mission weight fraction

Wn_W0 = (warmup * taxi * takeoff * climb * cruise * combat * loiter * descent * goaround * goaround * landing)       # 2 goaround allowance per the RFP requirements

# Fuel fractions

F_used = 1 - Wn_W0          # Used fuel fraction
F = 1.06 * F_used           # Total fuel fraction with reserves

# Initial Weight Guess and Convergence Loop

Wo = 75000                  # Initial guess for Takeoff Gross Weight in lbs
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


# Plot Convergence
import matplotlib.pyplot as plt
# Plot convergence
plt.figure(figsize=(8, 4))
plt.title("Weight Estimate Convergence")
plt.xlabel("Iteration")
plt.ylabel("W₀ (kg)")
plt.plot(Wo_history, label="W₀", linewidth=2)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Print results from weight estimation
print("Gross Takeoff Weight: " + str(round(Wo)) + " lbs")
print("Empty Weight: " + str(round(We)) + " lbs")
print("Regression's Empty Weight Fraction (We/W0): " + str(round(empty_weight_fraction ,3)))
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
