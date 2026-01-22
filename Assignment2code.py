# Stage 1 Inputs
Wcrew = 400
Wpayload = 6802

# Regression constants
A = 2.392
c = -.13
Kvs = 1.04  # Variable sweep wing constant (1.04 for variable, 1.0 for fixed)

# Mission segment weight fractions
warmup   = 0.99
taxi     = 0.99
takeoff  = 0.99
climb    = 0.96
cruise   = 0.812
combat   = 0.90
loiter   = 0.977
descent  = 0.99
landing  = 0.995

# Mission weight fraction
Wn_W0 = (warmup * taxi * takeoff * climb * cruise * combat * loiter * descent * landing)

# Fuel fractions
F_used = 1 - Wn_W0          # mission fuel
F = 1.06 * F_used           # mission + reserve

# Initial guess
Wo = 75000
previous_Wo = 0


for iteration in range(30):
    E = A * Wo**c * Kvs
    denominator = 1 - F - E
    
    Wo_new = (Wcrew + Wpayload) / denominator

    if abs(Wo_new - Wo) < 1:
        print(f"Converged at iteration {iteration+1}")
        break

    Wo = Wo_new

E = min(A * Wo**c * Kvs, 0.51)

# Final weights
We = E * Wo # Empty Weight
Wfuel_total = F * Wo # Total fuel weight
Wfuel_used = F_used * Wo # Used Fuel weight
Wfuel_reserved = Wfuel_total - Wfuel_used # Reserve Fuel weight
W_landing = Wn_W0 * Wo # Landing weight

# ============================
# Cost Estimating Relationships
# ============================

def engineering_cost(W_airframe, V_H, Q, F_cert, F_CF, F_comp, F_press, F_HyE, R_ENG, CPI):
    C_ENG = (0.083 * W_airframe**0.791 * V_H**1.521 * Q**0.183 * F_cert * F_CF * F_comp * F_press * F_HyE * R_ENG * CPI)
    return C_ENG

def tooling_cost(W_airframe, V_H, Q, Q_M, F_taper, F_CF, F_comp, F_press, F_HyE, R_TOOL, CPI):
    C_TOOL = (2.1036 * W_airframe**0.764 * V_H**0.899 * Q**0.178 * Q_M**0.066 * F_taper * F_CF * F_comp * F_press * F_HyE * R_TOOL * CPI)
    return C_TOOL

def manufacturing_cost(W_airframe, V_H, Q, F_cert, F_CF, F_comp, F_HyE, R_MFG, CPI):
    C_MFG = (20.2588 * W_airframe**0.74 * V_H**0.543 * Q**0.524 * F_cert * F_CF * F_comp * F_HyE * R_MFG * CPI)
    return C_MFG

def development_support_cost(W_airframe, V_H, Q_proto, F_cert, F_CF, F_comp, F_press, F_HyE, CPI):
    C_DEV = (0.06458 * W_airframe**0.873 * V_H**1.89 * Q_proto**0.346 * F_cert * F_CF * F_comp * F_press * F_HyE * CPI)
    return C_DEV

def flight_test_operations_cost(W_airframe, V_H, Q_proto, F_cert, F_HyE, CPI):
    C_FT = (0.009646 * W_airframe**1.16 * V_H**1.3718 * Q_proto**1.281 * F_cert * F_HyE * CPI)
    return C_FT

# ============================
# CER Parameters (5th Gen Fighter Estimates)
# ============================
W_airframe = We                 # Empty weight (lbs)
V_H = 1190                      # Max velocity (knots) - F-35 reference
Q = 500                        # Total production quantity
Q_M = 100                       # Manufacturing rate (per year)
Q_proto = 5                     # Prototype quantity

# Adjustment Factors (typical for 5th gen fighter)
F_cert = 1.3                    # Certification factor
F_CF = 0.95                     # Commonality factor
F_comp = 1.25                   # Complexity factor
F_press = 1.15                  # Pressure vessel factor
F_taper = 0.95                  # Taper factor
F_HyE = 1.2                     # High yield equipment factor

# Recurring cost indices
R_ENG = 1.0                     # Engine recurring
R_TOOL = 1.0                    # Tooling recurring
R_MFG = 1.0                     # Manufacturing recurring

# Cost Price Index (escalation factor, 1.0 = baseline)
CPI = 1.0

# Calculate costs
C_ENG = engineering_cost(W_airframe, V_H, Q, F_cert, F_CF, F_comp, F_press, F_HyE, R_ENG, CPI)
C_TOOL = tooling_cost(W_airframe, V_H, Q, Q_M, F_taper, F_CF, F_comp, F_press, F_HyE, R_TOOL, CPI)
C_MFG = manufacturing_cost(W_airframe, V_H, Q, F_cert, F_CF, F_comp, F_HyE, R_MFG, CPI)
C_DEV = development_support_cost(W_airframe, V_H, Q_proto, F_cert, F_CF, F_comp, F_press, F_HyE, CPI)
C_FT = flight_test_operations_cost(W_airframe, V_H, Q_proto, F_cert, F_HyE, CPI)
unit_cost = (C_TOOL + C_MFG) / Q
C_RDTE = (C_ENG + C_DEV + C_FT)
C_TOTAL = (C_RDTE + C_TOOL + C_MFG)


print(f"\nFinal Takeoff Gross Weight: {Wo:.2f} lbs")
print(f"Empty Weight Fraction: {E:.4f}")
print(f"Empty Weight: {We:.2f} lbs")
print(f"Landing Weight: {W_landing:.2f} lbs")
print(f"Total Fuel Weight: {Wfuel_total:.2f} lbs")
print(f"Used Fuel Weight: {Wfuel_used:.2f} lbs")
print(f"Reserve Fuel Weight: {Wfuel_reserved:.2f} lbs")

print(f"\n--- Cost Estimates ---")
print(f"Engineering Cost: ${C_ENG:,.0f}")
print(f"Development Support Cost: ${C_DEV:,.0f}")
print(f"Flight Test Cost: ${C_FT:,.0f}")
print(f"RDT&E Total: ${C_RDTE:,.0f}")
print(f"\nTooling Cost: ${C_TOOL:,.0f}")
print(f"Manufacturing Cost: ${C_MFG:,.0f}")
print(f"\nTotal Program Cost: ${C_TOTAL:,.0f}")
print(f"Unit Production Cost: ${unit_cost:,.0f}")
