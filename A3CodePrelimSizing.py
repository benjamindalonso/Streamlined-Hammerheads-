# This script calculates important preliminary sizing parameters and outputs a plot showing a fesable design region

import numpy as np
import matplotlib.pyplot as plt
import math

# Create wing loading entries for x axis of plot
Wing_Loading = np.linspace(0, 100, 500)  # in lb/ft^2

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

# Maneuvering Constraint
def Maneuvering_Constraint(TurnRate, g, Vturn, Cd0Turn, Wing_Loading, k, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction):
    n = math.sqrt((((TurnRate * Vturn)/g)**2) + 1)
    q = (0.5 * rhoTurn * (Vturn ** 2))
    Wing_LoadingTurn = Wing_Loading * TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction
    ThrustToWeightTurn = ((q * Cd0Turn) / (Wing_LoadingTurn)) + ((k / q) * Wing_LoadingTurn * (n ** 2))
    Maneuvering_Constraint = ThrustToWeightTurn * ((TakeoffFuelFraction * ClimbFuelFraction)/(ThrustReduction))
    return Maneuvering_Constraint

# Launch Constraint
# Insert Launch Constraint def Here

# Landing Constraint
# Insert Landing Constraint def Here

# Ceiling Constraint
# Insert Ceiling Constraint def Here

# Dash Constraint
# Insert Dash Constraint def Here

# Climb Constraint
# Insert Climb Constraint def Here


# PARAMETERS
# Add additional parameters to the bottom of the list as needed
Clmax = 1.5 # Maximum coefficient of lift (this will occur right at stall - max angle of attack)
rho = 0.0023769 # Air density at stall condition in slugs/ft^3
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




# CALCULATIONS
Cruise = Cruise_Constraint(rhoCruise, Vcruise, Cd0Cruise, k, Wing_Loading, TakeoffFuelFraction, ClimbFuelFraction,ThrustReduction)
Stall = Stall_Constraint(Clmax, rho, Vstall)
Maneuver = Maneuvering_Constraint(TurnRate, g, Vturn, Cd0Turn, Wing_Loading, k, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction)
#Launch = Launch_Constraint()  # Fill in parameters
#Landing = Landing_Constraint()  # Fill in parameters
#Ceiling = Ceiling_Constraint()  # Fill in parameters
#Dash = Dash_Constraint()  # Fill in parameters
#Climb = Climb_Constraint()  # Fill in parameters



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
plt.plot(Wing_Loading, Maneuver, color='green', linewidth=2.5, 
         label='Maneuvering Constraint')


# Add lines here to plot the constraint you added




# Formatting
plt.xlim(0, 50)                # Typical fighter range
plt.ylim(0, 1.2)                # T/W usually 0.8–1.2 for fighters
plt.xlabel('Wing Loading  W/S  (psf)', fontsize=12)
plt.ylabel('Thrust-to-Weight Ratio  T/W', fontsize=12)
plt.title('Preliminary Sizing Constraint Diagram\nNext-Gen Carrier-Based Strike Fighter', fontsize=14, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='upper left', fontsize=11, framealpha=0.9)

plt.tight_layout()
plt.show()