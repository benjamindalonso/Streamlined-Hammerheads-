"""
Aircraft Component Weight Estimation (Raymer Equations)
"""
import numpy as np

# ==========================
# INPUT PARAMETERS
# ==========================

# Create weight and position arrays
weights = np.array([])
weights_empty = np.array([])
X_positions = np.array([])

# General aircraft parameters
W_dg      = 58619.0     # Design gross weight (lbs)
N_z       = 8.0         # Ultimate load factor
M         = 0.8         # Cruise Mach number
N_c       = 1           # Number of crew
thrust    = 43000       # Maximum Thrust (lbf)
S_fw      = 30.915      # Firewall area (ft^2)
L_ec      = 22          # Length of engine wiring (ft)
N_s       = 3           # Number of flight control systems
N_u       = 6           # Number of hydraulic functions
K_mc      = 1.45        # mission completion mulitplier for electrical systems
R_kva     = 110         # Voltage rating for aircraft (kV)
L_a       = 27          # dist. of generators -> avionics -> cockpit (ft)
N_gen     = 3           # number of generators on engine
TSFC      = .889        # Thrust specific fuel consumption in lb/lbf/hr (Assumed value for a modern fighter engine at cruise)

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
S_vt      = 43.10771   # Vertical tail area (ft²) (for only one)
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
D_e       = 3.583       # Engine diameter (feet)
N_en      = 1           # Number of engines
K_vg      = 1.0         # Variable geometry intake factor (1.0 for fixed geometry)
K_d       = 1.0         # Duct type / design factor (1.0 for our case)

# Fuel system parameters
N_t = 5                 # number of fuel tanks
V_t = 420               # total fuel volume (ft^2)
V_i = 0                 # integral tank volume
V_p = V_t               # volume self sealing tanks (We're assuming none of the tanks are integrated into stucture, for now)

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

weights = np.append(weights,[W_Engine,W_ForwardTank,W_MainTank,W_DropTank,3*W_AIM9,W_AIM120,W_Avionics])
weights_empty = np.append(weights_empty,[W_Engine,0,0,0,0,0,W_Avionics])

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

X_firewall = 15 # between cockpit and forward fuel tank
X_engine_controls = 25.754 # assumed to be at front of engine
X_flight_controls = 24 # assumed to be near middle of fuselage
X_hydraulics = 24 # assumed to be near middle of fuselage
X_electrical = 24 # assumed to be near middle of fuselage
X_furnishings = 12.244 # CG of seat
X_AC = 24 # CG of airconditioning and anti ice (assumed to be near middle of fuselage)
X_handling_gear = 24 # CG assumed to be near middle of fuselage
X_fuel_system_and_tanks = 32.939 # Assumed to be at CG of fuel CG
X_pilot = 14.211

X_positions = np.append(X_positions,[X_Engine,X_ForwardTank,X_MainTank,X_DropTank,X_AIM9,X_120,X_Avionics])

# ==========================
# CALCULATIONS
# ==========================

import math

# Pilot and gear
W_pilot = 250
weights = np.append(weights,W_pilot)
weights_empty = np.append(weights_empty,W_pilot)
X_positions = np.append(X_positions,X_pilot)

# 15.1 Wing Weight
W_wing = (0.0103 * K_dw * K_vs * (W_dg * N_z)**0.5 * S_w**0.622 * AR_w**0.785 *
          (t_over_c_root)**(-0.4) * (1 + taper_ratio)**0.05 *
          (math.cos(math.radians(Lambda_w)))**(-1.0) * S_csw**0.04)
weights = np.append(weights,W_wing) # adds weight value to array
weights_empty = np.append(weights_empty,W_wing)
X_positions = np.append(X_positions,X_Wing)

# 15.2 Horizontal Tail Weight
W_horizontal_tail = (3.316 * (1 + F_w / B_h)**(-2.0) *
                     ((W_dg * N_z) / 1000)**0.260 * S_ht**0.806)
weights = np.append(weights,W_horizontal_tail)
weights_empty = np.append(weights_empty,W_horizontal_tail)
X_positions = np.append(X_positions,X_Ht)

# 15.3 Vertical Tail Weight
W_vertical_tail = 2*(0.452 * K_rht * (1 + H_t / H_v)**0.5 *
                   (W_dg * N_z)**0.488 * S_vt**0.718 * M**0.341 *
                   L_t**(-1.0) * (1 + S_r_over_S_vt)**0.348 * A_vt**0.223 *
                   (1 + taper_ratio)**0.25 * (math.cos(math.radians(Lambda_vt)))**(-0.323)) # 15.3 multiplied by two for the two vert tails
weights = np.append(weights,W_vertical_tail)
weights_empty = np.append(weights_empty,W_vertical_tail)
X_positions = np.append(X_positions,X_Vt)

# 15.4 Fuselage Weight
W_fuselage = (0.499 * K_dwf * W_dg**0.35 * N_z**0.25 * L_fuse**0.5 *
              D_fuse**0.849 * W**0.685)
weights = np.append(weights,W_fuselage)
weights_empty = np.append(weights_empty,W_fuselage)
X_positions = np.append(X_positions,X_Fuselage)

# 15.5 Main Landing Gear Weight
W_main_landing_gear = K_cb * K_tpg * (W_l * N_l)**0.25 * L_m**0.973
weights = np.append(weights,W_main_landing_gear)
weights_empty = np.append(weights_empty,W_main_landing_gear)
X_positions = np.append(X_positions,X_main_landing_gear)

# 15.6 Nose Landing Gear Weight
W_nose_landing_gear = (W_l * N_l)**0.290 * L_n**0.5 * N_nw**0.525
weights = np.append(weights,W_nose_landing_gear)
weights_empty = np.append(weights_empty,W_nose_landing_gear)
X_positions = np.append(X_positions,X_nose_landing_gear)

# 15.7 Engine Mounts
W_engine_mounts = 0.013*((thrust)**(0.579))*N_z
weights = np.append(weights,W_engine_mounts)
weights_empty = np.append(weights_empty,W_engine_mounts)
X_positions = np.append(X_positions,X_Engine) # Assume same as engine CG

# 15.8 Firewall
W_firewall = 1.135*S_fw
weights = np.append(weights,W_firewall)
X_positions = np.append(X_positions,X_firewall)
weights_empty = np.append(weights_empty,W_firewall)

# 15.9 Engine Section
W_engine_section = 0.01*(W_Engine**0.717)*N_z
weights = np.append(weights,W_engine_section)
weights_empty = np.append(weights_empty,W_engine_section)
X_positions = np.append(X_positions,X_Engine) # Assume same as engine CG

# 15.10 Air Induction System Weight (the intake)
W_air_induction = (13.29 * K_vg * L_d**0.643 * K_d**0.182 * 
                   N_en**1.498 * (L_s / L_d)**(-0.373) * D_e)
weights = np.append(weights,W_air_induction)
weights_empty = np.append(weights_empty,W_air_induction)
X_positions = np.append(X_positions,X_Intake)

# 15.13 Oil Cooling
W_oil_cooling = 37.28*N_en**1.023
weights = np.append(weights,W_oil_cooling)
weights_empty = np.append(weights_empty,W_oil_cooling)
X_positions = np.append(X_positions,X_Engine) # Assume same as engine CG

# 15.14 Engine Controls
W_engine_controls = 10.5*(N_en**1.008)*(L_ec**0.222)
weights = np.append(weights,W_engine_controls) 
weights_empty = np.append(weights_empty,W_engine_controls)
X_positions = np.append(X_positions,X_engine_controls) # Assume CG at front of engine

# 15.16 Fuel Systems and Tanks
W_fuel_system_and_tanks = 7.45*(V_t**0.47)*((1+V_i/V_t)**-0.095)*(1+V_p/V_t)*(N_t**0.06)*N_en*(((thrust*TSFC)/1000)**0.249)
weights = np.append(weights,W_fuel_system_and_tanks)
weights_empty = np.append(weights_empty,W_fuel_system_and_tanks)
X_positions = np.append(X_positions,X_fuel_system_and_tanks) # Assumed to be at CG of fuel stores

# 15.17 Flight Controls
W_flight_controls = 36.28*(M**0.003)*(S_csw**0.489)*(N_s**0.484)*(N_c**0.127)
weights = np.append(weights,W_flight_controls)
X_positions = np.append(X_positions,X_flight_controls) # Assume CG near middle of fuselage
weights_empty = np.append(weights_empty,W_flight_controls)

# 15.19 Hydraulics
W_hydraulics = 37.23*N_u**0.604
weights = np.append(weights,W_hydraulics)
X_positions = np.append(X_positions,X_hydraulics) # Assume CG near middle of fuselage
weights_empty = np.append(weights_empty,W_hydraulics)

# 15.20 Electrical
W_electrical = 172.2*K_mc*(R_kva**0.152)*(N_c**0.10)*(L_a**0.10)*(N_gen**0.091)
weights = np.append(weights,W_electrical)
X_positions = np.append(X_positions,X_electrical) # Assume CG near middle of fuselage
weights_empty = np.append(weights_empty,W_electrical)

# 15.22 Furnishings
W_furinishings = 217.6*N_c
weights = np.append(weights,W_furinishings)
X_positions = np.append(X_positions,X_furnishings)
weights_empty = np.append(weights_empty,W_furinishings)

# 15.23 Air Conditioning and Anti-Ice
W_AC = 201.6*((1400+200*N_c)/1000)**0.735
weights = np.append(weights,W_AC)
weights_empty = np.append(weights_empty,W_AC)
X_positions = np.append(X_positions,X_AC) # Assume CG near middle of fuselage

# 15.24 Handling Gear
W_handling_gear = (3.2*10**(-4))*W_dg
weights = np.append(weights,W_handling_gear)
weights_empty = np.append(weights_empty,W_handling_gear)
X_positions = np.append(X_positions,X_handling_gear)

X_Cg_Aircraft = np.sum(np.multiply(X_positions,weights))/np.sum(weights)
X_Cg_Aircraft_NoFuelorArms = np.sum(np.multiply(X_positions,weights_empty))/np.sum(weights_empty)

# Fuel Fraction Variables
TSFC = .889 # Thrust specific fuel consumption in lb/lbf/hr (Assumed value for a modern fighter engine at cruise) Also defined on line 29
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
print(f"Engine Mount Weight      = {W_engine_mounts:.2f} lbs")
print(f"Firewall Weight          = {W_firewall:.2f} lbs")
print(f"Engine Section Weight    = {W_engine_section:.2f} lbs")
print(f"Air Induction System     = {W_air_induction:.2f} lbs")   # Added
print(f"Oil Colling              = {W_oil_cooling:.2f} lbs")
print(f"Engine Controls          = {W_engine_controls:.2f} lbs")
print(f"Fuel System and Tanks    = {W_fuel_system_and_tanks:.2f} lbs")
print(f"Flight Controls          = {W_flight_controls:.2f} lbs")
print(f"Hydraulics               = {W_hydraulics:.2f} lbs")
print(f"Electrical               = {W_electrical:.2f} lbs")
print(f"Furnishings              = {W_furinishings:.2f} lbs")
print(f"AC and Anti-ice          = {W_AC:.2f} lbs")
print(f"Handling Gear            = {W_handling_gear:.2f} lbs")
print(f"X_Cg_Aircraft            = {X_Cg_Aircraft:.2f} ft")
print(f"X_Cg_Aircraft_NoFuelorArms = {X_Cg_Aircraft_NoFuelorArms:.2f} ft")



total_emptyish = np.sum(weights_empty)
total_weight = np.sum(weights)

print(f"\nSum of these components = {total_emptyish:.1f} lbs")
print(f"GTOW                    = {total_weight:.2f} lbs")
print (f"\nStartup & Taxi Fuel Fraction = {StartUp_and_Taxi:.4f}")
print (f"Takeoff Fuel Fraction = {Takeoff:.4f}")
print (f"Loiter Fuel Fraction = {Loiter:.4f}")