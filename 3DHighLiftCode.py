# Variables

S_Slats = 210.884 # Surface area of wing affected by slats
S_Flaps = 285.2 # Surface area of wing affected by flaps
S_ref = 600 # Surface area of main wing
Delta_Cl_Flap = .9 # Flap Clmax approximation
Delta_Cl_Slat = .3 # Slat Clmax approximation
Clmax_Clean = 1.8793 # Clmax for clean configuration
WingSweep = 38 # Wing sweep angle
FlapSweep = 82.1 # Flap sweep angle
SlatSweep = 46.1 # Slat Sweep angle
CosWingSweep = .788
CosFlapSweep = .137
CosSlatSweep = .693

# Equations

CL_Max_Clean = .9*(Clmax_Clean)*CosWingSweep
CL_Max_Flaps = .9*(Delta_Cl_Flap)*(S_Flaps/S_ref)*CosFlapSweep
CL_Max_Slats = .9*(Delta_Cl_Slat)*(S_Slats/S_ref)*CosSlatSweep

# Print

print("CLMaxClean =", round(CL_Max_Clean, 5), "Unitless")
print("CLMaxFlaps =", round(CL_Max_Flaps, 5), "Unitless")
print("CLMaxSlats =", round(CL_Max_Slats, 5), "Unitless")