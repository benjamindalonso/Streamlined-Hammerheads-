# This script calculates important preliminary sizing parameters and outputs a plot showing a fesable design region

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
    TW_turn = (q * Cd0 / Wing_LoadingTurn) + (k * (n**2) / q * Wing_LoadingTurn)
    
    # Correct back to Sea Level Static T/W (T0/W0)
    # (T/W)_turn * (Weight_at_maneuver / W0) * (1 / Thrust_Reduction)
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
def Landing_Constraint(GTOW,landWeight,maxLandSpeed,CLmaxLand,density):
    S_land = 2*landWeight/(density*((maxLandSpeed/1.15)**2)*CLmaxLand)
    landWingLoading = GTOW/S_land
    return landWingLoading

# Ceiling Constraint
def Ceiling_Constraint(Cd0Cruise,k):
    absolute_ceiling = np.ones_like(Wing_Loading) * (2 * np.sqrt(k * Cd0Cruise))
    return absolute_ceiling

# Insert Ceiling Constraint def Here

# Dash Constraint
# Insert Dash Constraint def Here

# Climb Constraint
# Insert Climb Constraint def Here


# PARAMETERS
# Add additional parameters to the bottom of the list as needed
Clmax = 1.5 # Maximum coefficient of lift (this will occur right at stall - max angle of attack)
rho = 0.0023769 # Air density at stall condition in slugs/ft^3
rhoTropicalDay = 0.00219 # Air density at sea level on a tropical day in slugs/ft^3 (for launch constraint)
Vstall = 135 # Airspeed at stall 
rhoCruise = 0.0007382 # Air density at cruise altitude in slugs/ft^3
Vcruise = 550 # Cruise velocity in knots
Cd0Cruise = 0.01 # Zero lift drag coefficient at cruise
k = 1.2 # Stall speed factor 
TakeoffFuelFraction = 0.99 # Fuel fraction from assignment2code
ClimbFuelFraction = 0.96 # Fuel fraction from assignment2code
ThrustReduction = 0.8 # Thrust reduction factor at cruise (due to altitude and speed) Get this from engine data
TurnRate = 0.1745 # in Rad/s based on rfp preference of 10 deg/s
g = 32.174 # ft/s^2 (gravitational acceleration)
Vturn = 500 # Velocity during the maneuver in feet per second
Cd0Turn = 0.01 # Zero lift drag coefficient during the turn (I just used the openVSP value again)
rhoTurn = 0.001267 # Slugs per ft^3 (Density at 20,000 ft - maneuvering altitude per rfp)
MidMissionFuelFraction = 0.906 # Fuel fraction half way through cruise portion of mission (This is based on the rfp requirements for maneuvering being done at "mid mission weight")
Vend = 135 # Catipult end speed in knots with a 67,000 GTOW and a 210 CSV setting on the catipult 
Vwod = 0 # Wind speed over the deck in knots (Assumed 0 for worst case scenario)
Vthrust = 10 # Velocity added by engine thrust during catipult launch (Assumed to be 10 knots per Raymer page 136)
ClmaxTakeOff = 1.7 # Clmax at takeoff per slide 11 of preliminary sizing part 2


# CALCULATIONS
Cruise = Cruise_Constraint(rhoCruise, Vcruise, Cd0Cruise, k, Wing_Loading, TakeoffFuelFraction, ClimbFuelFraction,ThrustReduction)
Stall = Stall_Constraint(Clmax, rho, Vstall)
#Launch = Launch_Constraint()  # Fill in parameters
Launch = Launch_Constraint(rhoTropicalDay, Vend, Vwod, Vthrust, ClmaxTakeOff) 
Landing = Landing_Constraint(67822,51010,202.6,1.5,23.77*10**(-4))
Ceiling = Ceiling_Constraint(Cd0Cruise,k)
#Dash = Dash_Constraint()  # Fill in parameters
#Climb = Climb_Constraint()  # Fill in parameters
Maneuver_Sustained = Sustained_Turn_Constraint(TurnRate, g, Vturn, Cd0Turn, k, Wing_Loading, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction)
Maneuver_Instant = Instantaneous_Turn_Constraint(Clmax, TurnRate, Vturn, g, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction)



# PLOTTING
plt.figure(figsize=(10, 7))
plt.style.use('ggplot')  # or 'seaborn' / 'default' for cleaner look

# Plot stall constraint (vertical line)
plt.axvline(x=Stall, color='red', linestyle='--', linewidth=2.5, 
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
plt.axvline(x=Launch, color='yellow', linewidth=2.5, 
            label=f'Launch Constraint (W/S ≤ {Launch:.1f} psf)')

# Plot landing constraint
plt.axvline(x=Landing,color='pink',linestyle='-.', linewidth=2.5,
            label='Landing Constraint (W/S to left permissible)')

# Plot ceiling consstraint
plt.plot(Wing_Loading,Ceiling, color='black', linewidth=2.5, label= 'ceilingconstraint')




# Formatting
plt.xlim(0, 100)                # Typical fighter range
plt.ylim(0, 1.2)                # T/W usually 0.8–1.2 for fighters
plt.xlabel('Wing Loading  W/S  (psf)', fontsize=12)
plt.ylabel('Thrust-to-Weight Ratio  T/W', fontsize=12)
plt.title('Preliminary Sizing Constraint Diagram\nNext-Gen Carrier-Based Strike Fighter', fontsize=14, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='upper left', fontsize=11, framealpha=0.9)

plt.tight_layout()
plt.show()