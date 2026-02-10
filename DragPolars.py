# This code estimates the drag polars (cl vs cd) for 5 different flight configurations
import numpy as np
import matplotlib.pyplot as plt


# DEFINITIONS

# Estimate Swet
def Swet_Estimator(Wto, c, d):
    Swet = (10 ** c) * (Wto ** d) 
    return Swet

# Estiamte Equivilant Parasite Area (f)
def f_Estimator(Cf, Swet):
    f = Cf * Swet
    return f

# Estimate clean (no flaps, gear up) zero lift coefficient of drag (slide 52 of sizing lecture part 1)
def Cd0_Estimator(Swet, Sref, a, b):
    f = (a * Swet ** b)
    Cd0 = f / Sref
    return Cd0, f

# Function for calculating drag polar
def DragPolar(Cl, Cd0, AR, e):
    Cd = Cd0 + (Cl ** 2) / (np.pi * AR * e)
    return Cd




# PARAMETERS
 
Sref = 740 # ft^2 (Based on openVSP area of main wing + cannard)
a = -2.3979 # Slide 52 of sizing lecture part 1
b = 1 # Slide 52 of sizing lecture part 1
Cf = 0.0040 # Table 12.3 (Navy Fighter)
Wto = 67822 # lbf From Assignment 2 Weight Estimate
c = -0.1289 # Table 3.5 (Fighter)
d = 0.7506  # Table 3.5 (Fighter)
AR = 2.5 # Based on openVSP model
e_Clean = 0.8 # We need to refine this estimation based on historical data (0.8 is just a guess)
Cd0_Clean = 0.01 # Based on openVSP drag analysis 
e_TakeOffFlaps = e_Clean * 0.8 # Based on slide 56 of sizing lecture part 1
Cd0_TakeOffFlaps = Cd0_Clean + 0.02 # Based on slide 56 of sizing lecture part 1
e_LandingFlaps = e_Clean * 0.75 # Based on slide 56 of sizing lecture part 1
Cd0_LandingFlaps = Cd0_Clean + 0.075 # Based on
e_TakeOffFlapsGearDown = e_TakeOffFlaps # Gear doesn't effect wing efficiency 
Cd0_TakeOffFlapsGearDown = Cd0_TakeOffFlaps + 0.025 # Based on slide 56 of sizing lecture part 1
e_LandingFlapsGearDown = e_LandingFlaps # Gear doesn't effect wing efficiency
Cd0_LandingFlapsGearDown = Cd0_LandingFlaps + 0.025 # Based on slide 56 of sizing lecture part 1
Swet = Swet_Estimator(Wto, c, d)
f = f_Estimator(Cf, Swet)
Cl = np.linspace(-3, 3, 100) 



# CALCULATIONS

# Calculate Cd for clean configuration
Cd_Clean = DragPolar(Cl, Cd0_Clean, AR, e_Clean)

# Calculate Cd for Take Off Flaps Down Gear Up
Cd_TakeOffFlaps = DragPolar(Cl, Cd0_TakeOffFlaps, AR, e_TakeOffFlaps)

# Calculate Cd for Landing Flaps Down Gear Up
Cd_LandingFlaps = DragPolar(Cl, Cd0_LandingFlaps, AR, e_LandingFlaps)

# Calculate Cd for Take Off Flaps Down Gear Down
Cd_TakeOffFlapsGearDown = DragPolar(Cl, Cd0_TakeOffFlapsGearDown, AR, e_TakeOffFlapsGearDown)

# Calculate Cd for Landing Flaps Down Gear Down
Cd_LandingFlapsGearDown = DragPolar(Cl, Cd0_LandingFlapsGearDown, AR, e_LandingFlapsGearDown)



# PLOTTING
 
plt.figure(figsize=(10, 7))

# Clean
plt.plot(Cd_Clean, Cl, color='royalblue', linewidth=2, label='Clean')

# Takeoff flaps, gear up
plt.plot(Cd_TakeOffFlaps, Cl, color='limegreen', linewidth=2, label='Takeoff flaps, gear up')

# Landing flaps, gear up
plt.plot(Cd_LandingFlaps, Cl, color='darkorange', linewidth=2, label='Landing flaps, gear up')

# Takeoff flaps + gear down
plt.plot(Cd_TakeOffFlapsGearDown, Cl, color='crimson', linewidth=2.0, label='Takeoff flaps + gear down')

# Landing flaps + gear down
plt.plot(Cd_LandingFlapsGearDown, Cl, color='purple', linewidth=2.0, label='Landing flaps + gear down')

# Formatting
plt.xlabel('CD  (drag coefficient)', fontsize=12)
plt.ylabel('CL  (lift coefficient)', fontsize=12)
plt.title('Drag Polar – Different Configurations', fontsize=13)

plt.grid(True, alpha=0.3, linestyle='--')

# Reasonable axis limits (adjust if needed)
plt.xlim(0.0, 0.9)      # most configs should fit here
plt.ylim(-3, 3)      # allows seeing negative CL if desired

plt.legend(loc='upper left', fontsize=10, framealpha=0.95)

plt.tight_layout()
plt.show()

print('Swet:', Swet)