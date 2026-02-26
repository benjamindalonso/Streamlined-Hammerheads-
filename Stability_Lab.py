import math
import numpy as np

# =========================================================
# GEOMETRY INPUTS
# =========================================================

S = 543.04                # Wing area (ft^2)
b = 35.83                 # Wing span (ft)
MAC = 17.07               # Mean Aerodynamic Chord (ft)

Sh = 68.93                # Horizontal tail area (ft^2)
Htail_arm = 15.29         # Distance from CG to horizontal tail (ft)

Sv = 86.22                # Vertical tail area (ft^2)

Fuselage_length = 50      # Fuselage length (ft)
Max_Fuselage_width = 11.82

CG_Location = 29.19       # From nose (ft)
Wing_25_MAC_along_CRL = 4.27  # From nose (ft)

# =========================================================
# TAIL VOLUME COEFFICIENTS
# =========================================================

Cht = (Sh * Htail_arm) / (S * MAC)
Cvt = (Sv * Htail_arm) / (S * b)

print("\n================ STABILITY ANALYSIS ================")
print("Horizontal Tail Volume Coefficient (Cht):", round(Cht,4))
print("Vertical Tail Volume Coefficient (Cvt):", round(Cvt,4))

# =========================================================
# ASPECT RATIOS
# =========================================================

AR_w = b**2 / S
b_h = 20  # <-- Replace with actual horizontal tail span (ft)
AR_h = b_h**2 / Sh

print("\nWing Aspect Ratio (AR_w):", round(AR_w,3))
print("Horizontal Tail Aspect Ratio (AR_h):", round(AR_h,3))

# =========================================================
# LIFT CURVE SLOPES
# =========================================================

eta = 0.97
M = 0.85
Lambda = 0

a_w_exact = (2 * math.pi * AR_w) / (
    2 + math.sqrt(((AR_w/eta)**2) * (1 + math.tan(Lambda)**2 - M**2) + 4)
)

a_h_exact = (2 * math.pi * AR_h) / (
    2 + math.sqrt(((AR_h/eta)**2) * (1 + math.tan(Lambda)**2 - M**2) + 4)
)

print("\nWing Lift Curve Slope (per rad):", round(a_w_exact,3))
print("Clean Tail Lift Curve Slope (per rad):", round(a_h_exact,3))

# =========================================================
# DOWNWASH
# =========================================================

downwash_exact = (2 / (math.pi * AR_w)) * a_w_exact

print("\nDownwash Gradient (dε/dα):", round(downwash_exact,3))

# =========================================================
# FUSELAGE DESTABILIZING TERM
# =========================================================

K_f = 0.115

dCm_dCL_fus = (
    (K_f * Max_Fuselage_width * Fuselage_length)
    / (S * MAC)
) * (1 / a_w_exact)

print("Fuselage Contribution (dCm_fuselage/dCL):", round(dCm_dCL_fus,4))

# =========================================================
# EFFECTIVE TAIL LIFT CURVE SLOPE
# =========================================================

a_h_effective = a_h_exact * (1 - downwash_exact)

print("Effective Tail Lift Curve Slope (per rad):", round(a_h_effective,3))

# =========================================================
# NEUTRAL POINT
# =========================================================

hn = 0.25 + (a_h_effective / a_w_exact) * Cht - dCm_dCL_fus

print("\nNeutral Point (fraction of MAC):", round(hn,4))
print("Neutral Point (% MAC):", round(hn*100,2))

# =========================================================
# CG LOCATION
# =========================================================

h_cg = ((CG_Location - Wing_25_MAC_along_CRL) / MAC) + 0.25

print("\nCG Location (fraction of MAC):", round(h_cg,4))
print("CG Location (% MAC):", round(h_cg*100,2))

# =========================================================
# STATIC MARGIN
# =========================================================

static_margin_percent = (hn - h_cg) * 100

print("\nStatic Margin (% MAC):", round(static_margin_percent,2))
print("======================================================")