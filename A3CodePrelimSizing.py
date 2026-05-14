import numpy as np
import matplotlib.pyplot as plt
import math

# Create wing loading entries for x axis of plot
Wing_Loading = np.linspace(0.5, 100, 500)  # in lb/ft^2

# Stall Constraint (This is a verticle line on the wingloading/thrustloading plot
def Stall_Constraint(Clmax, rho, Vstall):
    Stall_Constraint = 0.5 * rho * (Vstall**2) * Clmax
    return Stall_Constraint

# Cruise Constraint
def Cruise_Constraint(rhoCruise, Vcruise, Cd0Cruise, k, Wing_Loading, TakeoffFuelFraction, ClimbFuelFraction,ThrustReduction):
    qcr = (0.5 * rhoCruise * (Vcruise ** 2))
    Wing_LoadingCruise = Wing_Loading * TakeoffFuelFraction * ClimbFuelFraction
    ThrustToWeightCruise = ((qcr * Cd0Cruise) / (Wing_LoadingCruise)) + ((k / qcr) * Wing_LoadingCruise)
    Cruise_Constraint = ThrustToWeightCruise * ((TakeoffFuelFraction * ClimbFuelFraction)/(ThrustReduction))
    return Cruise_Constraint

# Maneuvering Constraint (sustained and instant)
def Sustained_Turn_Constraint(TurnRate, g, Vturn, Cd0, k, Wing_Loading, rho, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction):
    n = math.sqrt(((TurnRate * Vturn) / g)**2 + 1)
    q = 0.5 * rho * (Vturn**2)
    Wing_LoadingTurn = Wing_Loading * TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction 
    
    # Required T/W at maneuver alt
    TW_turn = (q * Cd0 / Wing_LoadingTurn) + (K * (n**2) * Wing_LoadingTurn / q )
    
    # Correct back to Sea Level Static T/W (T0/W0)
    TW_SLS = TW_turn * (TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction / ThrustReduction)
    return TW_SLS

def Instantaneous_Turn_Constraint(Clmax, TurnRate, Vturn, g, rho, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction):
    n = math.sqrt(((TurnRate * Vturn) / g)**2 + 1)
    q = 0.5 * rho * (Vturn**2)
    WS_limit = (q * Clmax) / (n * TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction)
    return WS_limit

# Launch Constraint
def Launch_Constraint(rhoTropicalDay, Vend, Vwod, Vthrust, ClmaxTakeOff):
    Launch_Constraint = (0.5) * rhoTropicalDay * ((Vend + Vwod + Vthrust)**2) * (ClmaxTakeOff) / 1.21
    return Launch_Constraint

# Landing Constraint finds wing loading for landing weights lbf,speeds ft/s, density slug/ft^3
def Landing_Constraint(GTOW, landWeight, maxLandSpeed, CLmaxLand, density):
    S_land = 2*landWeight/(density*((maxLandSpeed/1.15)**2)*CLmaxLand)
    landWingLoading = GTOW/S_land
    return landWingLoading

# Ceiling Constraint
def Ceiling_Constraint(Cd0Cruise, K):
    absolute_ceiling = np.ones_like(Wing_Loading) * (2 * np.sqrt(K * Cd0Cruise))
    return absolute_ceiling

# Dash Constraint
def Dash_Constraint(rhoDash, MachDash, aDash, CD0Dash, K, Wing_Loading, ThrustReduction):
    Vdash = MachDash * aDash             
    qdash = 0.5 * rhoDash * Vdash**2
    TW_dash_at_alt = (qdash * CD0Dash) / Wing_Loading + (K / qdash) * Wing_Loading
    TW_SLS = TW_dash_at_alt / ThrustReduction
    return TW_SLS

# Climb Constraint
def Climb_Constraint(Ks, K, Climb_Cd0, Clmax, Climb_Gradient):
    Climb_Intial = ((Ks**2*Climb_Cd0)/(Clmax))+((K*((Clmax)/Ks**2))+((Climb_Gradient)))
    Climb_Constraint = Climb_Intial*((1/.8)*(1/.99))
    return Climb_Constraint


# PARAMETERS
Clmax = 1.5
rho = 0.0023769
rhoTropicalDay = 0.00219
Vstall = 135
rhoCruise = 0.0007382
Vcruise = 550
Cd0Cruise = 0.00696
Ks = 1.8
TakeoffFuelFraction = 0.99
ClimbFuelFraction = 0.96
ThrustReduction = 0.8
TurnRate = 0.1745
g = 32.174
Vturn = 500
Cd0Turn = 0.00696
rhoTurn = 0.001267
MidMissionFuelFraction = 0.906
ROC = 200
V_horizontal = 500 * 1.68781
V_horizontal_min = V_horizontal * 60
Climb_Gradient = ROC / V_horizontal_min
Climb_Cd0 = 0.01696
e = 0.8
AR = 2.5
K = 1/(math.pi*e*AR)
rhoDash = 0.000889
aDash   = 994.0
MachDash = 2
CD0Dash  = Cd0Cruise
Vend = 135
Vwod = 0
Vthrust = 10
ClmaxTakeOff = 1.
GTOW = 67822
LandingWeight = 51010
maxLandSpeed = 202.6
ClmaxLand = 1.5

# Design Point Parameters
design_weight = 54748.1  # lbf
wing_area = 600          # ft^2
total_thrust = 43000     # lbf

design_WS = design_weight / wing_area
design_TW = total_thrust / design_weight


# CALCULATIONS
Cruise = Cruise_Constraint(rhoCruise, Vcruise, Cd0Cruise, K, Wing_Loading, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction)
Stall = Stall_Constraint(Clmax, rho, Vstall)
Launch = Launch_Constraint(rhoTropicalDay, Vend, Vwod, Vthrust, ClmaxTakeOff)
Landing = Landing_Constraint(GTOW, LandingWeight, maxLandSpeed, ClmaxLand, rhoTropicalDay)
Ceiling = Ceiling_Constraint(Cd0Cruise, K)
Dash = Dash_Constraint(rhoDash, MachDash, aDash, CD0Dash, K, Wing_Loading, ThrustReduction)
Climb = Climb_Constraint(Ks, K, Climb_Cd0, Clmax, Climb_Gradient)
Maneuver_Sustained = Sustained_Turn_Constraint(TurnRate, g, Vturn, Cd0Turn, K, Wing_Loading, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction)
Maneuver_Instant = Instantaneous_Turn_Constraint(Clmax, TurnRate, Vturn, g, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction)


# PLOTTING
plt.figure(figsize=(10, 7))
plt.style.use('default')

# Plot stall constraint (vertical line)
plt.axvline(x=Stall, color='red', linewidth=2.5,
            label=f'Stall Constraint (W/S ≤ {Stall:.1f} psf)')

# Plot cruise constraint
plt.plot(Wing_Loading, Cruise, color='blue', linewidth=2.5,
         label='Cruise Constraint')

# Plot maneuvering constraint
plt.plot(Wing_Loading, Maneuver_Sustained, color='green', linewidth=2.5,
         label='Maneuvering Constraint (sustained)')
plt.axvline(x=Maneuver_Instant, color='purple', linestyle='-.', linewidth=2.5,
            label=f'Instant Maneuvering Constraint (W/S ≤ {Maneuver_Instant:.1f})')

# Plot launch constraint (vertical line)
plt.axvline(x=Launch, color='black', linestyle='-.', linewidth=2.5,
            label=f'Launch Constraint (W/S ≤ {Launch:.1f} psf)')

# Plot landing constraint
plt.axvline(x=Landing, color='red', linestyle='-.', linewidth=2.5,
            label='Landing Constraint (W/S to left permissible)')

# Plot climb constraint
plt.axhline(y=Climb, color='brown', linestyle='-.', linewidth=2.5,
            label='Climb Constraint (T/W ≥ {:.2f})'.format(Climb))

# Plot ceiling constraint
plt.plot(Wing_Loading, Ceiling, color='black', linewidth=2.5,
         label='Ceiling Constraint')

# Plot dash constraint
plt.plot(Wing_Loading, Dash, color='orange', linewidth=1,
         label='Dash Constraint')

# Plot design point
plt.scatter(design_WS, design_TW,
            color='r', marker='*', s=300,
            label=f'Design Point (W/S={design_WS:.1f} psf, T/W={design_TW:.2f})')

# Formatting
plt.xlim(0, 100)
plt.ylim(0, 1.2)
plt.xlabel('Wing Loading  W/S  (psf)', fontsize=12)
plt.ylabel('Thrust-to-Weight Ratio  T/W', fontsize=12)
plt.title('Preliminary Sizing Constraint Diagram\nNext-Gen Carrier-Based Strike Fighter', fontsize=14, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='lower right', fontsize=11, framealpha=0.9)

plt.tight_layout()
plt.show()