# Parameters
FuelDensity = 6.7 # lbm/gal
AdditionalFuelBuffer = 0.15 # 15% additional fuel for safety
FuelWeight = 19000 # lbf (from A4 code)


# Print Results
print(f"Required Fuel Mass Without Safety Buffer: {RequiredFuel:.2f} lbm")
print(f"Required Fuel Volume With Safety Buffer: {RequiredFuelVolumeAdjusted:.2f} gal")