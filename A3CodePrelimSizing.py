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
    TW_turn = (q * Cd0 / Wing_LoadingTurn) + (K * (n**2) / q * Wing_LoadingTurn)
    
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

def Landing_Constraint(GTOW,landWeight,maxLandSpeed,CLmaxLand,density):
    S_land = 2*landWeight/(density*((maxLandSpeed/1.15)**2)*CLmaxLand)
    landWingLoading = GTOW/S_land
    return landWingLoading

# Ceiling Constraint

def Ceiling_Constraint(Cd0Cruise,K):
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
    Climb_Constraint = Climb_Intial*((1/.8)*(1/.99)) # Adjusting for fuel fractions and thrust reduction
    return Climb_Constraint



# PARAMETERS

# Add additional parameters to the bottom of the list as needed
Clmax = 1.5 # Maximum coefficient of lift (this will occur right at stall - max angle of attack)
rho = 0.0023769 # Air density at stall condition in slugs/ft^3
rhoTropicalDay = 0.00219 # Air density at sea level on a tropical day in slugs/ft^3 (for launch constraint)
Vstall = 135 # Airspeed at stall in knots  
rhoCruise = 0.0007382 # Air density at cruise altitude in slugs/ft^3
Vcruise = 550 # Cruise velocity in knots
Cd0Cruise = 0.00696 # Zero lift drag coefficient at cruise
Ks = 1.8 # Stall speed factor 
TakeoffFuelFraction = 0.99 # Fuel fraction from assignment2code
ClimbFuelFraction = 0.96 # Fuel fraction from assignment2code
ThrustReduction = 0.8 # Thrust reduction factor at cruise (due to altitude and speed) Get this from engine data
TurnRate = 0.1745 # in Rad/s based on rfp preference of 10 deg/s
g = 32.174 # ft/s^2 (gravitational acceleration)
Vturn = 500 # Velocity during the maneuver in feet per second
Cd0Turn = 0.00696 # Zero lift drag coefficient during the turn (I just used the openVSP value again)
rhoTurn = 0.001267 # Slugs per ft^3 (Density at 20,000 ft - maneuvering altitude per rfp)
MidMissionFuelFraction = 0.906 # Fuel fraction half way through cruise portion of mission (This is based on the rfp requirements for maneuvering being done at "mid mission weight")
ROC = 200 # Rate of climb ft/min
V_horizontal = 500 * 1.68781  # Climb airspeed knots to ft/s
V_horizontal_min = V_horizontal * 60  # Climb airspeed in ft/min
Climb_Gradient = ROC / V_horizontal_min # Climb gradient
Climb_Cd0 = 0.01696 # Climb drag coefficient (I just used the openVSP value again)
e = 0.8 # Oswald efficiency factor (typical value for fighters)
AR = 2.5 # Aspect ratio (typical value for fighters)
K = 1/(math.pi*e*AR) # Induced drag factor
rhoDash = 0.000889 # Density at 30,000 ft in slugs/ft^3
aDash   = 994.0 # Speed of sound at 30k ft in ft/s
MachDash = 2 #Dash speed in Mach
CD0Dash  = Cd0Cruise # Assuming zero lift drag coefficient at dash is the same as cruise
Vend = 135 # Catipult end speed in knots with a 67,000 GTOW and a 210 CSV setting on the catipult 
Vwod = 0 # Wind speed over the deck in knots (Assumed 0 for worst case scenario)
Vthrust = 10 # Velocity added by engine thrust during catipult launch (Assumed to be 10 knots per Raymer page 136)
ClmaxTakeOff = 1.7 # Clmax at takeoff per slide 11 of preliminary sizing part 2
GTOW = 67822 # Gross takeoff weight from assignment 2
LandingWeight = 51010 # Landing weight based on assignment 2 fuel fractions 
maxLandSpeed = 202.6 # Max landing speed in feet per second 
ClmaxLand = 1.5 # Maximum CL during landing phase 


# CALCULATIONS
Cruise = Cruise_Constraint(rhoCruise, Vcruise, Cd0Cruise, K, Wing_Loading, TakeoffFuelFraction, ClimbFuelFraction,ThrustReduction)
Stall = Stall_Constraint(Clmax, rho, Vstall)
#Launch = Launch_Constraint()  # Fill in parameters
Launch = Launch_Constraint(rhoTropicalDay, Vend, Vwod, Vthrust, ClmaxTakeOff) 
Landing = Landing_Constraint(GTOW,LandingWeight,maxLandSpeed,ClmaxLand,rhoTropicalDay)
Ceiling = Ceiling_Constraint(Cd0Cruise, K)
Dash = Dash_Constraint(rhoDash, MachDash, aDash, CD0Dash, K, Wing_Loading, ThrustReduction)
# Baseline climb = 45,000 ft/min
Climb = Climb_Constraint(Ks, K, Climb_Cd0, Clmax, Climb_Gradient)
#Climb = Climb_Constraint()  # Fill in parameters
Maneuver_Sustained = Sustained_Turn_Constraint(TurnRate, g, Vturn, Cd0Turn, K, Wing_Loading, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction)
Maneuver_Instant = Instantaneous_Turn_Constraint(Clmax, TurnRate, Vturn, g, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction)



# PLOTTING
plt.figure(figsize=(10, 7))
plt.style.use('ggplot')  # or 'seaborn' / 'default' for cleaner look

# Plot stall constraint (vertical line)
plt.axvline(x=Stall, color='red',  linewidth=2.5, 
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
plt.axvline(x=Launch, color='black', linewidth=2.5, 
            label=f'Launch Constraint (W/S ≤ {Launch:.1f} psf)')

# Plot landing constraint
plt.axvline(x=Landing,color='red',linestyle='-.', linewidth=2.5,
            label='Landing Constraint (W/S to left permissible)')

# Plot climb constraint
plt.axhline(y=Climb, color='brown', linestyle='-.', linewidth=2.5,
            label='Climb Constraint (T/W ≥ {:.2f})'.format(Climb))
# Plot ceiling consstraint
plt.plot(Wing_Loading,Ceiling, color='black', linewidth=2.5, label= 'Ceiling Constraint')

# Plot dash constraint
plt.plot(Wing_Loading, Dash, color='orange', linewidth=2.5, label='Dash Constraint')

# Formatting

plt.xlim(0, 100)                # Typical fighter range
plt.ylim(0, 1.2)                # T/W usually 0.8–1.2 for fighters
plt.xlabel('Wing Loading  W/S  (psf)', fontsize=12)
plt.ylabel('Thrust-to-Weight Ratio  T/W', fontsize=12)
plt.title('Preliminary Sizing Constraint Diagram\nNext-Gen Carrier-Based Strike Fighter', fontsize=14, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='lower right', fontsize=11, framealpha=0.9)

plt.tight_layout()

plt.fill_between(Wing_Loading,1.2,Dash,where=(Wing_Loading<=Launch),color='lightsteelblue', label='Feasible Region')
plt.show()