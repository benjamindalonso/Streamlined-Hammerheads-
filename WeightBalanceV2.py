"""
Aircraft Component Weight Estimation (Raymer Equations)
"""

# ==========================
# INPUT PARAMETERS
# ==========================

# General aircraft parameters
W_dg      = 58619.0     # Design gross weight (lbs)
N_z       = 8.0         # Ultimate load factor
M         = 0.8         # Cruise Mach number

# Wing parameters
S_w       = 600.0       # Wing reference area (ft²)
AR_w      = 2.028       # Wing aspect ratio
t_over_c_root = 0.06    # Root thickness-to-chord ratio 
Lambda_w   = 38.0       # Wing sweep angle (degrees) 
taper_ratio = 0.24      # Wing taper ratio (λ)
K_dw       = 1.0        # Wing weight factor for Delta Wing (Page 610 Raymer)
K_vs       = 1.0        # Variable sweep factor (1.0 for fixed sweep Page 610 Raymer))
quarterMAC = 21.54008   # Distance from nose to quarter MAC on main wing in feet
MAC = 19.35           # Mean Aerodynamic Chord (ft) 
Span = 34.8          # Wing span (ft)   

# Horizontal tail parameters
S_ht      = 77.27030    # Horizontal tail area (ft²)
F_w       = 6.361        # Fuselage width at Canard intersection (ft)
B_h       = 22.972    # Horizontal tail span (ft)

# Vertical tail parameters
quarterMACTail = 38.7313 # Distance from nose to quarter MAC on tail in feet
S_vt      = 86.21542    # Vertical tail area (ft²)
A_vt      = 2.88328     # Vertical tail aspect ratio
H_t       = 0.001       # Horizontal tail height above fuselage reference line (ft) (Page 610 Raymers 1.0 for T tail, 0 for conventional tail)
H_v       = 0.001      # Vertical tail height above fuselage reference line (ft) (Page 610 Raymers 1.0 for T tail, 0 for conventional tail)
L_t       = quarterMACTail - quarterMAC  # Wing quarter MAC to vertical tail quarter MAC in ft
S_r_over_S_vt = 0.3   # Rudder area / vertical tail area (Sr/Svt)
Lambda_vt = 36.35714    # Vertical tail sweep angle (degrees)
K_rht     = 1.0         # Vertical tail weight correction factor (usually 1.0)


# Fuselage parameters
L_fuse    = 47.78        # Fuselage length (ft) - you had 49.5 earlier
D_fuse    = 11.41728    # Fuselage maximum diameter or equivalent width (ft) 
K_dwf     = 1.0         # Fuselage weight correction factor (usually 1.0)
W = 13.5                # Fuselage structural width (widest point) in feet

# Landing gear parameters
W_l       = 58619.0     # Landing weight 
N_l       = 4.5         # Ultimate landing load factor (Page 611 Raymer Ngear * 1.5)
L_m       = 108         # Extended main landing gear length (inches) - strut length will be adjusted based on gear sizing
L_n       = 6.0         # Extended nose landing gear length (inches) - strut length will be adjusted based on gear sizing
N_nw      = 2           # Number of nose wheels
K_cb      = 1.0         # Landing gear weight correction (1.0 for conventional)
K_tpg     = 1.0         # Usually 1 (page 610 Raymer)

# Air Induction System (Intake) parameters - NEW
L_d       = 8.16429     # Length of diffuser / air induction system (ft)
L_s       = 2           # Length of straight intake duct (ft)
D_e       = 3.583        # Engine diameter (feet)
N_en      = 1           # Number of engines
K_vg      = 1.0         # Variable geometry intake factor (1.0 for fixed geometry)
K_d       = 1.0         # Duct type / design factor (1.0 for our case)

# Control Surfaces
S_csw     =   ((0.3 * MAC) * (0.175 * Span)) + ((0.4 * MAC) * (0.175 * Span))      # Control surface area includes flaps (wing mounted) 

# Other Weights 
W_Engine = 6422 # lbs
W_ForwardTank = 3007
W_MainTank = 15035.86   
W_DropTank = 3007
W_AIM9 = 380
W_AIM120 = 696
W_Avionics = 2500

# Locations of component CGs from nose in feet
X_Intake = 24.776
X_Engine = 40.576
X_Vt = 44.849
X_Ht = 18.476
X_Wing = 31.505
X_ForwardTank = 23.16
X_MainTank = 37.500
X_DropTank = 34.683
X_AIM9 = 37.163
X_120 = 38.67
X_Avionics = 7.276
X_Fuselage = 27.104
X_nose_landing_gear = 15.853  # Estimated location of nose landing gear CG from nose (ft)
X_main_landing_gear = 40.164  # Estimated location of main landing gear



# ==========================
# CALCULATIONS
# ==========================

import math

# 15.1 Wing Weight
W_wing = (0.0103 * K_dw * K_vs * (W_dg * N_z)**0.5 * S_w**0.622 * AR_w**0.785 *
          (t_over_c_root)**(-0.4) * (1 + taper_ratio)**0.05 *
          (math.cos(math.radians(Lambda_w)))**(-1.0) * S_csw**0.04)

# 15.2 Horizontal Tail Weight
W_horizontal_tail = (3.316 * (1 + F_w / B_h)**(-2.0) *
                     ((W_dg * N_z) / 1000)**0.260 * S_ht**0.806)

# 15.3 Vertical Tail Weight
W_vertical_tail = (0.452 * K_rht * (1 + H_t / H_v)**0.5 *
                   (W_dg * N_z)**0.488 * S_vt**0.718 * M**0.341 *
                   L_t**(-1.0) * (1 + S_r_over_S_vt)**0.348 * A_vt**0.223 *
                   (1 + taper_ratio)**0.25 * (math.cos(math.radians(Lambda_vt)))**(-0.323))

# 15.4 Fuselage Weight
W_fuselage = (0.499 * K_dwf * W_dg**0.35 * N_z**0.25 * L_fuse**0.5 *
              D_fuse**0.849 * W**0.685)

# 15.5 Main Landing Gear Weight
W_main_landing_gear = K_cb * K_tpg * (W_l * N_l)**0.25 * L_m**0.973

# 15.6 Nose Landing Gear Weight
W_nose_landing_gear = (W_l * N_l)**0.290 * L_n**0.5 * N_nw**0.525

# 15.10 Air Induction System Weight (the intake)
W_air_induction = (13.29 * K_vg * L_d**0.643 * K_d**0.182 * 
                   N_en**1.498 * (L_s / L_d)**(-0.373) * D_e)

X_Cg_Aircraft = ((W_fuselage*X_Fuselage)+(W_vertical_tail*X_Vt)+(W_horizontal_tail*X_Ht)+(W_wing*X_Wing)+(W_air_induction*X_Intake)+(W_Engine*X_Engine)+(W_ForwardTank*X_ForwardTank)+(W_MainTank*X_MainTank)+(W_DropTank*X_DropTank)+(W_AIM9*X_AIM9)+(W_AIM120*X_120)+(W_AIM120*X_120)+(W_AIM120*X_120)+(W_Avionics*X_Avionics)+(W_nose_landing_gear * X_nose_landing_gear)+(W_main_landing_gear * X_main_landing_gear))/((W_fuselage)+(W_vertical_tail)+(W_horizontal_tail)+(W_wing)+(W_air_induction)+(W_Engine)+(W_ForwardTank)+(W_MainTank)+(W_DropTank)+(W_AIM9)+(W_AIM120)+(W_AIM120)+(W_AIM120)+(W_Avionics)+(W_nose_landing_gear)+(W_main_landing_gear))
X_Cg_Aircraft_NoFuelorArms = ((W_fuselage*X_Fuselage)+(W_vertical_tail*X_Vt)+(W_horizontal_tail*X_Ht)+(W_wing*X_Wing)+(W_air_induction*X_Intake)+(W_Engine*X_Engine)+(W_Avionics*X_Avionics))/((W_fuselage)+(W_vertical_tail)+(W_horizontal_tail)+(W_wing)+(W_air_induction)+(W_Engine)+(W_Avionics))

# Fuel Fraction Variables
TSFC = .889 # Thrust specific fuel consumption in lb/lbf/hr (Assumed value for a modern fighter engine at cruise)
TSFC_SeaLevel = .3 # Thrust specific fuel consumption at sea level in lbm/lbf/hr (Assumed value for a modern fighter engine at sea level)
TSFC_Loiter = .7 # Thrust specific fuel consumption during loiter in lb/lbf/hr (Assumed value for a modern fighter engine during loiter)
TSFC_Cruise = .8 # Thrust specific fuel consumption during cruise in lb/lbf/hr (Assumed value for a modern fighter engine during cruise)
T_Max = 43000 # Maximum thrust in pounds (Assumed value for a modern fighter engine)
T_Idle = .05 * T_Max # Idle thrust in pounds (Assumed to be 5% of max thrust)
t_Idle = 15 # Time spent Idle in minutes (Assumed value for startup and taxi)
t_Idle_Hours = t_Idle / 60 # Time spent Idle in hours
t_Takeoff = 1 # Time spent in takeoff in minutes (Assumed value for takeoff roll and initial climb)
t_Takeoff_Hours = t_Takeoff / 60 # Time spent in takeoff in hours
Max_L_D = 13 # Maximum lift-to-drag ratio for the aircraft (Assumed value for a modern fighter)
E = 20 # Loiter time in minutes (Assumed value for a modern fighter)
t_Loiter_Hours = E / 60 # Time spent in loiter in hours
# Fuel Fractions

StartUp_and_Taxi = 1 - (t_Idle_Hours * TSFC_SeaLevel * (T_Idle / W_dg)) # Fuel fraction for startup and taxi
Takeoff = 1 - (t_Takeoff_Hours * TSFC * (T_Max / W_dg))
Loiter = math.exp((-t_Loiter_Hours * TSFC_Loiter) / Max_L_D) # Fuel fraction for loiter


# ==========================
# OUTPUT
# ==========================

print(f"Wing Weight              = {W_wing:.2f} lbs")
print(f"Horizontal Tail Weight   = {W_horizontal_tail:.2f} lbs")
print(f"Vertical Tail Weight     = {W_vertical_tail:.2f} lbs")
print(f"Fuselage Weight          = {W_fuselage:.2f} lbs")
print(f"Main Landing Gear Weight = {W_main_landing_gear:.2f} lbs")
print(f"Nose Landing Gear Weight = {W_nose_landing_gear:.2f} lbs")
print(f"Air Induction System     = {W_air_induction:.2f} lbs")   # Added
print(f"X_Cg_Aircraft = {X_Cg_Aircraft:.2f} ft")
print(f"X_Cg_Aircraft_NoFuelorArms = {X_Cg_Aircraft_NoFuelorArms:.2f} ft")



total_emptyish = (W_wing + W_horizontal_tail + W_vertical_tail + W_fuselage +
                  W_main_landing_gear + W_nose_landing_gear + W_air_induction)

print(f"\nSum of these components = {total_emptyish:.1f} lbs")
print (f"Startup & Taxi Fuel Fraction = {StartUp_and_Taxi:.4f}")
print (f"Takeoff Fuel Fraction = {Takeoff:.4f}")
print (f"Loiter Fuel Fraction = {Loiter:.4f}")