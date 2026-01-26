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
empty_weight_fraction_percent = empty_weight_fraction * 100  # Empty weight fraction in percent

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





# Cost Estimating Relationships


def engineering_cost(W_airframe, V_H, Q, E_cert, E_CF, E_comp, E_press, E_HyE, R_ENG, CPI): # Engineering Cost CER
    C_ENG = (0.083 * W_airframe**0.791 * V_H**1.521 * Q**0.183 * E_cert * E_CF * E_comp * E_press * E_HyE * R_ENG * CPI) # Engineering Cost Calculation
    return C_ENG # Return Engineering Cost

def tooling_cost(W_airframe, V_H, Q, Q_M, T_taper, T_CF, T_comp, T_press, T_HyE, R_TOOL, CPI): # Tooling Cost CER
    C_TOOL = (2.1036 * W_airframe**0.764 * V_H**0.899 * Q**0.178 * Q_M**0.066 * T_taper * T_CF * T_comp * T_press * T_HyE * R_TOOL * CPI) # Tooling Cost Calculation
    return C_TOOL # Return Tooling Cost

def manufacturing_cost(W_airframe, V_H, Q, M_cert, M_CF, M_comp, M_HyE, R_MFG, CPI): # Manufacturing Cost CER
    C_MFG = (20.2588 * W_airframe**0.74 * V_H**0.543 * Q**0.524 * M_cert * M_CF * M_comp * M_HyE * R_MFG * CPI) # Manufacturing Cost Calculation
    return C_MFG # Return Manufacturing Cost

def development_support_cost(W_airframe, V_H, Q_proto, D_cert, D_CF, D_comp, D_press, D_HyE, CPI): # Development Support Cost CER
    C_DEV = (0.06458 * W_airframe**0.873 * V_H**1.89 * Q_proto**0.346 * D_cert * D_CF * D_comp * D_press * D_HyE * CPI) # Development Support Cost Calculation
    return C_DEV # Return Development Support Cost

def flight_test_operations_cost(W_airframe, V_H, Q_proto, F_cert, F_HyE, CPI): # Flight Test Operations Cost CER
    C_FT = (0.009646 * W_airframe**1.16 * V_H**1.3718 * Q_proto**1.281 * F_cert * F_HyE * CPI) # Flight Test Operations Cost Calculation
    return C_FT # Return Flight Test Operations Cost


# CER Parameters (5th Gen Fighter Estimates)

W_airframe = .45*(We)           # Estimate between 40-50% of empty weight
V_H = 1303.46                   # Max velocity (knots) - F-22 reference
Q = 500                         # Total production quantity over 5 years
Q_M = 100                       # Manufacturing rate (per year)
Q_proto = 5                     # Prototype quantity

# Adjustment Factors for Engineering Cost

E_cert = .67                  # Certification factor for LSA Certification
E_CF = 1.03                   # Commonality factor for complex flap system
E_comp = 2                    # Complexity factor for 100% Composites
E_press = 1.03                # Pressure vessel factor for pressurized cabin
E_HyE = 1.47                  # High yield equipment factor for hybrid electric propulsion

# Adjustment Factors for Tooling Cost   

T_CF = 1.02                  # Commonality factor for complex flap system
T_comp = 2                   # Complexity factor for 100% Composites
T_press = 1.01               # Pressure vessel factor for pressurized cabin
T_Taper = .95                # Factor for untapered wings
T_HyE = 1.1                  # High yield equipment factor for hybrid electric propulsion

# Adjustment Factors for Manufacturing Cost

M_cert = .75                 # Certification factor for LSA Certification
M_CF = 1.01                  # Commonality factor for complex flap system
M_comp = 1.25                # Complexity factor for 100% Composites
M_HyE = 1.1                  # High yield equipment factor for hybrid electric propulsion

# Adjustment Factors for Development Support Cost

D_cert = .5                  # Certification factor for LSA Certification
D_CF = 1.01                  # Commonality factor for complex flap system
D_comp = 1.5                 # Complexity factor for 100% Composites
D_press = 1.03               # Pressure vessel factor for pressurized cabin
D_HyE = 1.05                 # High yield equipment factor for hybrid electric propulsion

# Adjustment Factors for FLight Test Operations Cost

F_cert = .5                  # Certification factor for LSA Certification
F_HyE = 1.5                  # High yield equipment factor for hybrid electric propulsion


# Recurring cost indices

R_ENG = (2.576*2026)-5058      # Engine recurring
R_TOOL = (2.883*2026)-5666     # Tooling recurring
R_MFG = (2.316*2026)-4552      # Manufacturing recurring

# Cost Price Index relative to base year (2012)

CPI = 1.41                     # Cost Price Index for 2026 relative to 2012

# Cost Calculations

C_ENG = engineering_cost(W_airframe, V_H, Q, F_cert, E_CF, E_comp, E_press, E_HyE, R_ENG, CPI) # Engineering Cost
C_TOOL = tooling_cost(W_airframe, V_H, Q, Q_M,T_Taper, T_CF, T_comp, T_press, T_HyE, R_TOOL, CPI) # Tooling Cost
C_MFG = manufacturing_cost(W_airframe, V_H, Q, M_cert, M_CF, M_comp, M_HyE, R_MFG, CPI) # Manufacturing Cost
C_DEV = development_support_cost(W_airframe, V_H, Q_proto, D_cert, D_CF, D_comp, D_press, D_HyE, CPI) # Development Support Cost
C_FT = flight_test_operations_cost(W_airframe, V_H, Q_proto, F_cert, F_HyE, CPI) # Flight Test Operations Cost
unit_cost = (C_TOOL + C_MFG) / Q        # Unit Production Cost
C_RDTE = (C_ENG + C_DEV + C_FT)         # RDT&E Total Cost
C_TOTAL = (C_RDTE + C_TOOL + C_MFG)     # Total Program Cost

# Weight Estimates Output
print(f"\n--- Weight Estimates ---") # Print Weight Estimates Header
print("Gross Takeoff Weight: " + str(round(Wo)) + " lbs")
print(f"Empty Weight: {We:.2f} lbs") # Print Empty Weight
print("Regression's Empty Weight Fraction (We/W0): " + str(round(empty_weight_fraction ,3)))
print(f"Regression's Empty Weight Percentage: {empty_weight_fraction:.2f}%") # Print Regression's Empty Weight Percentage
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



# Cost Estimates Output

print(f"\n--- Cost Estimates ---") # Print Cost Estimates Header
print(f"Engineering Cost: ${C_ENG:,.0f}") # Print Engineering Cost
print(f"Development Support Cost: ${C_DEV:,.0f}") # Print Development Support Cost
print(f"Flight Test Cost: ${C_FT:,.0f}") # Print Flight Test Cost
print(f"RDT&E Total: ${C_RDTE:,.0f}") # Print RDT&E Total Cost
print(f"\nTooling Cost: ${C_TOOL:,.0f}") # Print Tooling Cost
print(f"Manufacturing Cost: ${C_MFG:,.0f}") # Print Manufacturing Cost
print(f"\nTotal Program Cost: ${C_TOTAL:,.0f}") # Print Total Program Cost
print(f"Unit Production Cost: ${unit_cost:,.0f}") # Pr