# Parameters
Range = 2000 # nm
AvgSpeed = 550 # kts
SFC = 0.62 # lbm/hr/lbf
Thrust = 28000 # lbf
FuelDensity = 6.7 # lbm/gal

# Calculate Mission Time
MissionTime = Range / AvgSpeed # hr

# Calculate Required Fuel
RequiredFuel = SFC * MissionTime * Thrust # lbm (or lbf if at 1g and sea level)

# Calculate Required Fuel Volume
RequiredFuelVolume = RequiredFuel / FuelDensity # gal