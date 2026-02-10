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

# Landing Constraint finds wing loading for landing weights lbf,speeds ft/s, density slug/ft^3
def landing_Constraint(GTOW,landWeight,maxLandSpeed,CLmaxLand,density):
    S_land = 2*landWeight/(density*((maxLandSpeed/1.15)**2)*CLmaxLand)
    landWingLoading = GTOW/S_land
    return landWingLoading

# Ceiling Constraint
# Insert Ceiling Constraint def Here

# Dash Constraint
# Insert Dash Constraint def Here

# Climb Constraint

def Climb_Constraint(Ks, K, Climb_Cd0, Clmax, Climb_Gradient):
    Climb_Intial = ((Ks**2*Climb_Cd0)/(Clmax))+((K*((Clmax)/Ks**2))+((Climb_Gradient)))
    Climb_Constraint = Climb_Intial*((1/.8)*(1/.99)) # Adjusting for fuel fractions and thrust reduction
    return Climb_Constraint


# PARAMETERS
# Add additional parameters to the bottom of the list as needed
Clmax = 1.5 # Maximum coefficient of lift (this will occur right at stall - max angle of attack)
rho = 0.0023769 # Air density at stall condition in slugs/ft^3
Vstall = 135 # Airspeed at stall 
rhoCruise = 0.0007382 # Air density at cruise altitude in slugs/ft^3
Vcruise = 550 # Cruise velocity in knots
Cd0Cruise = 0.01 # Zero lift drag coefficient at cruise
Ks = 1.8 # Stall speed factor 
TakeoffFuelFraction = 0.99 # Fuel fraction from assignment2code
ClimbFuelFraction = 0.96 # Fuel fraction from assignment2code
ThrustReduction = 0.8 # Thrust reduction factor at cruise (due to altitude and speed) Get this from engine data
TurnRate = 0.1745 # in Rad/s based on rfp preference of 10 deg/s
g = 32.174 # ft/s^2 (gravitational acceleration)
Vturn = 500 # Velocity during the maneuver in feet per second
Cd0Turn = 0.01 # Zero lift drag coefficient during the turn (I just used the openVSP value again)
rhoTurn = 0.001267 # Slugs per ft^3 (Density at 20,000 ft - maneuvering altitude per rfp)
MidMissionFuelFraction = 0.906 # Fuel fraction half way through cruise portion of mission (This is based on the rfp requirements for maneuvering being done at "mid mission weight")
ROC = 200       # ft/min
V_horizontal = 500 * 1.68781  # knots to ft/s
V_horizontal_min = V_horizontal * 60  # ft/min
Climb_Gradient = ROC / V_horizontal_min
Climb_Cd0 = 0.068 # Climb drag coefficient (I just used the openVSP value again)
e = 0.8 # Oswald efficiency factor (typical value for fighters)
AR = 2.5 # Aspect ratio (typical value for fighters)
K = 1/(math.pi*e*AR) # Induced drag factor




# CALCULATIONS
Cruise = Cruise_Constraint(rhoCruise, Vcruise, Cd0Cruise, K, Wing_Loading, TakeoffFuelFraction, ClimbFuelFraction,ThrustReduction)
Stall = Stall_Constraint(Clmax, rho, Vstall)
Maneuver = Maneuvering_Constraint(TurnRate, g, Vturn, Cd0Turn, Wing_Loading, K, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction)
#Launch = Launch_Constraint()  # Fill in parameters
Landing = landing_Constraint(67822,51010,202.6,1.5,23.77*10**(-4))
#Ceiling = Ceiling_Constraint()  # Fill in parameters
#Dash = Dash_Constraint()  # Fill in parameters
# Baseline climb = 45,000 ft/min
Climb = Climb_Constraint(Ks, K, Climb_Cd0, Clmax, Climb_Gradient)



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

# Plot landing constraint
plt.axvline(x=Landing,color='red',linestyle='-.', linewidth=2.5,
            label='Landing Constraint (W/S to left permissible)')

# Plot climb constraint
plt.axhline(y=Climb, color='purple', linestyle='-.', linewidth=2.5,
            label='Climb Constraint (T/W ≥ {:.2f})'.format(Climb))

# Formatting
plt.xlim(0, 80)                # Typical fighter range
plt.ylim(0, 1.2)                # T/W usually 0.8–1.2 for fighters
plt.xlabel('Wing Loading  W/S  (psf)', fontsize=12)
plt.ylabel('Thrust-to-Weight Ratio  T/W', fontsize=12)
plt.title('Preliminary Sizing Constraint Diagram\nNext-Gen Carrier-Based Strike Fighter', fontsize=14, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='upper left', fontsize=11, framealpha=0.9)

plt.tight_layout()
plt.show()

print(Climb)
print(Climb_Gradient)