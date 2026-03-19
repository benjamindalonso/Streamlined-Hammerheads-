# Parameters
FuelDensity = 6.7 # lbf/gal
AdditionalFuelBuffer = 0.15 # 15% additional fuel for safety
FuelWeight = 13640 # lbf (from A4 code)

# Calculate Fuel Volume and Add Safety Buffer
RequiredFuelVolume = FuelWeight / FuelDensity
RequiredFuelVolumeAdjusted = RequiredFuelVolume * (1 + AdditionalFuelBuffer)


# Print Results
print(f"Required Fuel Volume With Safety Buffer: {RequiredFuelVolumeAdjusted:.2f} gal")