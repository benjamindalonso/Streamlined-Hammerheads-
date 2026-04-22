# Variables

S_Slats = 365.175 # Surface area of wing affected by slats
S_Flaps = 372.519 # Surface area of wing affected by flaps
S_ref = 600 # Surface area of main wing
Delta_Cl_Flap = .9  # Flap Clmax approximation
Delta_Cl_Slat = .3  # Slat Clmax approximation
Clmax_Clean = 1.8793 # Clmax for clean configuration
WingSweep = 38 # Wing sweep angle
FlapSweep = 82.1 # Flap sweep angle
SlatSweep = 46.1 # Slat Sweep angle
CosWingSweep = .788
CosFlapSweep = .994
CosSlatSweep = .723

# Equations 
CL_Max_Clean  = .9*(Clmax_Clean)*CosWingSweep

    
DeltaCL_Max_Flaps_Takeoff = .9*(Delta_Cl_Flap*.7)*(S_Flaps/S_ref)*CosFlapSweep
DeltaCL_Max_Slats_Takeoff = .9*(Delta_Cl_Slat*.7)*(S_Slats/S_ref)*CosSlatSweep


DeltaCL_Max_Flaps_Landing = .9*(Delta_Cl_Flap)*(S_Flaps/S_ref)*CosFlapSweep
DeltaCL_Max_Slats_Landing = .9*(Delta_Cl_Slat)*(S_Slats/S_ref)*CosSlatSweep

CL_Max_Takeoff = CL_Max_Clean+DeltaCL_Max_Flaps_Takeoff+DeltaCL_Max_Slats_Takeoff
CL_Max_Landing = CL_Max_Clean+DeltaCL_Max_Flaps_Landing+DeltaCL_Max_Slats_Landing
# Print

print("CLMaxTakeoff =", round(CL_Max_Takeoff, 5), "Unitless")
print("CLMaxLanding =", round(CL_Max_Landing, 5), "Unitless")
