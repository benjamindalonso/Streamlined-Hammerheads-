import numpy as np
import matplotlib.pyplot as plt
import math


#                        AIRCRAFT SIZING INPUTS


# Design Point
W0_design          = 54000.0      # Takeoff Gross Weight guess [lbf] (will be updated in loops)
S_design           = 600.0        # Wing area [ft²]
T_design           = 43000.0      # Total static thrust [lbf]

# Geometry & Configuration 
AR                 = 2.227        # Wing aspect ratio
e                  = 0.5398       # Oswald efficiency factor
K                  = 1 / (math.pi * e * AR)   # Induced drag factor

horiz_tail_area    = 77.7         # ft²
vert_tail_area     = 86.21        # ft²
wet_area_fuse      = 407.0        # fuselage wetted area [ft²]
number_of_engines  = 1

WeightCorrectionFactor = 0.87     # Scales empty weight to match known data point

# Aerodynamics 
Clean_Cd0          = 0.00766      # Zero-lift drag coefficient (clean)
Climb_Cd0          = Clean_Cd0 + 0.015   # Takeoff flaps

# Clmax values
Clmax_clean        = 1.6          # Clean configuration (for some constraints)
Clmax_takeoff      = 2.6
Clmax_landing      = 2.6
Clmax_stall        = 2.6
Ks                 = 1.8          # Stall speed safety factor for takeoff climb

# Atmosphere / Environment 
rho_SL             = 0.0023769    # Sea level standard [slug/ft³]
rho_tropical       = 0.00219      # Tropical day sea level [slug/ft³]
rho_20k            = 0.001267     # 20,000 ft [slug/ft³]
rho_30k            = 0.000891     # 30,000 ft [slug/ft³]
rho_cruise         = rho_30k

g                  = 32.174       # ft/s² (standard)

# Mission & Performance Requirements
# Speeds
Vstall_kts         = 120
Vstall_fps         = Vstall_kts * 1.68781

Vend_kts           = 135          # Catapult end speed
Vwod_kts           = 0            # Wind over deck
Vthrust_kts        = 10           # Thrust contribution during launch
Vend_fps           = Vend_kts * 1.68781
Vwod_fps           = Vwod_kts * 1.68781
Vthrust_fps        = Vthrust_kts * 1.68781

V_cruise_kts       = 500
V_cruise_fps       = V_cruise_kts * 1.68781

V_dash_kts         = 942          # Mach ~1.6 at 30k ft
V_dash_fps         = V_dash_kts * 1.68781

V_turn_kts         = 485
V_turn_fps         = V_turn_kts * 1.68781

V_landing_max_kts  = 140          # Max landing speed (approach)
V_landing_max_fps  = V_landing_max_kts * 1.68781

# Performance constraints
ROC_fpm            = 200          # Rate of climb [ft/min]
TurnRate_min       = 0.139626     # 8 deg/s [rad/s]
TurnRate_desired   = 0.174533     # 10 deg/s [rad/s]
n_struct_limit     = 8.0          # Structural g-limit for instantaneous turn

#  Fuel Fractions (Mission)
Wf_W0_fraction     = None         # Will be calculated in weight loop if used
TakeoffFuelFraction = 0.989
ClimbFuelFraction   = 0.96
MidMissionFuelFraction = 0.79     # At maneuver
CruiseFuelFraction  = TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction

#  Engine / Thrust Lapse 
ThrustReduction_cruise = 0.8      # T / T_SL at cruise altitude & speed
ThrustReduction_maneuver = 0.8
ThrustReduction_dash   = 0.8
Ct                 = 0.889        # TSFC [lb/hr/lbf]

#  Weights & Payload 
WeightCrew         = 250          # 1 pilot
WeightPayload      = (1000*4) + (190*2)    # MK83 JDAMs + AIM-9X 
WeightFuel = 21000   # Total Fuel Weight Based on A1 Analysis

#  Range / Endurance 
Range_nm           = 2000
Endurance_hr       = 0.333
MaxLD              = 13           # Maximum L/D (used for fuel fraction)

#  Constraint-Specific Parameters 
# Climb
Climb_Gradient     = (ROC_fpm) / (V_cruise_kts * 1.68781 * 60)   # sin(gamma) approx

# Cruise / Dash
q_cruise           = 0.5 * rho_30k * V_cruise_fps**2
q_dash             = 0.5 * rho_30k * V_dash_fps**2

# Maneuver
n_sustained        = None         # Calculated inside function

# Landing
LandingWeight = W0_design - (0.5 * WeightPayload) - (0.75 * WeightFuel) # Landing weight shall be with 25% fuel and 50% payload per RFP

# =============================================================================
#                        CONSTRAINT FUNCTIONS
# =============================================================================

def stall_ws_limit(rho=rho_tropical, Vstall_fps=Vstall_fps, Clmax=Clmax_stall):
    """Maximum allowable wing loading (vertical line) (To the left is feasible)"""
    return 0.5 * rho * Vstall_fps**2 * Clmax

def launch_ws_limit(rho=rho_tropical, Vend_fps=Vend_fps, Vwod_fps=Vwod_fps,
                    Vthrust_fps=Vthrust_fps, Clmax=Clmax_takeoff):
    """Maximum allowable wing loading from catapult launch (to the left is feasible)"""
    return 0.5 * rho * (Vend_fps + Vwod_fps + Vthrust_fps)**2 * Clmax / 1.21

def Landing_Constraint(W0_design = W0_design, LandingWeight = LandingWeight, V_landing_max_fps = V_landing_max_fps, Clmax_landing = Clmax_landing, rho_tropical = rho_tropical):
    S_land = 2*LandingWeight/(rho_tropical*((V_landing_max_fps/1.15)**2)*Clmax_landing)
    landWingLoading = W0_design/S_land
    return landWingLoading

def climb_tw_limit(Ks=Ks, Climb_Cd0=Climb_Cd0, Clmax=Clmax_takeoff,
                   K=K, Climb_Gradient=Climb_Gradient):
    """Minimum T/W (horizontal line)"""
    term1 = (Ks**2 * Climb_Cd0) / Clmax
    term2 = K * (Clmax / Ks**2)
    return (term1 + term2 + Climb_Gradient) * (1/0.8) * (1/0.99)

def cruise_tw(Wing_Loading, rho=rho_cruise, V_fps=V_cruise_fps, Cd0=Clean_Cd0,
              K=K, ff=CruiseFuelFraction, thrust_lapse=ThrustReduction_cruise):
    q = 0.5 * rho * V_fps**2
    WS_cruise = Wing_Loading * ff
    TW_alt = (q * Cd0 / WS_cruise) + (K * WS_cruise / q)
    return TW_alt * (ff / thrust_lapse)

def sustained_turn_tw(Wing_Loading, TurnRate=TurnRate_min, rho=rho_20k,
                      V_fps=V_turn_fps, Cd0=Clean_Cd0, K=K,
                      ff_mid=MidMissionFuelFraction,
                      ff_to=TakeoffFuelFraction, ff_climb=ClimbFuelFraction,
                      thrust_lapse=ThrustReduction_maneuver):
    n = math.sqrt(((TurnRate * V_fps) / g)**2 + 1)
    q = 0.5 * rho * V_fps**2
    ff = ff_to * ff_climb * ff_mid
    WS_turn = Wing_Loading * ff
    
    TW_turn = (q * Cd0 / WS_turn) + (K * (n**2) * WS_turn / q)
    return TW_turn * (ff / thrust_lapse)

def dash_tw(Wing_Loading, rho=rho_30k, V_fps=V_dash_fps, Cd0=Clean_Cd0,
            K=K, thrust_lapse=ThrustReduction_dash):
    q = 0.5 * rho * V_fps**2
    TW_alt = (q * Cd0 / Wing_Loading) + (K * Wing_Loading / q)
    return TW_alt / thrust_lapse   

# Instantaneous turn (max W/S limit)
def instant_turn_ws_limit(Clmax=Clmax_clean,                    # CLmax during maneuver (adjust as needed)
                          n_struct=n_struct_limit,                 # Structural limit +8g
                          rho=rho_20k,
                          V_fps=V_turn_fps):
    q = 0.5 * rho * V_fps**2
    
    # Wing loading at maneuver weight
    WS_mid = (q * Clmax) / n_struct
    
    # Convert to takeoff wing loading
    ff = TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction
    WS_takeoff = WS_mid / ff
    
    return WS_takeoff




# =============================================================================
#                        CALCULATIONS
# =============================================================================

# Create evaluation grid for Wing Loading (W/S)
Wing_Loading = np.linspace(5, 200, 800)   # psf

# Constraint Calculations
Stall_WS          = stall_ws_limit()
Launch_WS         = launch_ws_limit()
Climb_TW          = climb_tw_limit()
Cruise_TW         = cruise_tw(Wing_Loading)
Maneuver_Sustained_TW = sustained_turn_tw(Wing_Loading)
Maneuver_Instant_WS   = instant_turn_ws_limit()
Dash_TW           = dash_tw(Wing_Loading)
Landing_WS        = Landing_Constraint()



# Design Point (using your chosen values)
design_WS = W0_design / S_design
design_TW = T_design / W0_design


# =============================================================================
#                        PLOTTING - T/W vs W/S
# =============================================================================

plt.figure(figsize=(12, 8))
plt.style.use('default')

# Vertical constraints (maximum W/S)
plt.axvline(x=Stall_WS, color='red', linewidth=2.5, linestyle='--',
            label=f'Stall (W/S ≤ {Stall_WS:.1f} psf)')

plt.axvline(x=Launch_WS, color='black', linewidth=2.5, linestyle='-.',
            label=f'Launch (W/S ≤ {Launch_WS:.1f} psf)')

plt.axvline(x=Maneuver_Instant_WS, color='purple', linewidth=2.5, linestyle='-.',
            label=f'Instant Turn (+8g) (W/S ≤ {Maneuver_Instant_WS:.1f} psf)')

plt.axvline(x=Landing_WS, color='pink', linewidth=2.5, linestyle='--',
            label=f'Landing (W/S ≤ {Landing_WS:.1f} psf)')

# Curve constraints
plt.axhline(y=Climb_TW, color='brown', linewidth=2.5, linestyle='-.',
            label=f'Climb (T/W ≥ {Climb_TW:.3f})')

plt.plot(Wing_Loading, Cruise_TW, color='blue', linewidth=2.5, label='Cruise')
plt.plot(Wing_Loading, Maneuver_Sustained_TW, color='green', linewidth=2.5,
         label='Sustained Turn (8°/s)')
plt.plot(Wing_Loading, Dash_TW, color='orange', linewidth=2.5, label='Dash')



# =============================================================================
#                    COMPARABLE AIRCRAFT OVERLAY
# =============================================================================

aircraft_data = [
    ("F-22 Raptor",          840,  60000,  84000,  'k'),
    ("F-35 Lightning II",    460,  70000,  43000,  'm'),
    ("Rafale M",             492,  47000,  33720,  'g'),
    ("F/A-18E Super Hornet", 500,  66000,  44000,  'b'),
]

for name, S, TOGW, Total_T, color in aircraft_data:
    WS = TOGW / S
    TW = Total_T / TOGW
    plt.scatter(WS, TW, color=color, marker='o', s=80, edgecolors='black', 
                linewidth=1, zorder=6, label=name)

# Our Design Point on top
plt.scatter(design_WS, design_TW, color='r', marker='*', s=450, zorder=7,
            label=f'Our Design (W/S={design_WS:.1f}, T/W={design_TW:.3f})')

# ====================== SHADED FEASIBLE REGION ======================
x = Wing_Loading
y_lower = np.maximum(Maneuver_Sustained_TW, Dash_TW)
mask = x <= Maneuver_Instant_WS

plt.fill_between(x[mask], y_lower[mask], 1.6,
                 color='lightgreen', alpha=0.35,
                 label='Feasible Region')

# Formatting
plt.xlim(0, 200)
plt.ylim(0, 1.6)
plt.xlabel('Wing Loading  W/S  (lb/ft²)', fontsize=13)
plt.ylabel('Thrust-to-Weight Ratio  T/W', fontsize=13)
plt.title('Preliminary Sizing Constraint Diagram\n'
          'Next-Gen Carrier-Based Strike Fighter', fontsize=15, fontweight='bold')

plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='upper right', fontsize=10, framealpha=0.95)
plt.tight_layout()
plt.show()