# Reference Chord Length
reference_chord_length = 7  # feet

# Airspeed in feet per second
airspeed = 250  # ft/s

# Density of air at sea level in slugs/ft^3
density = 0.002243  # slugs/ft^3

# Dynamic viscosity of air at sea level in slugs/(ft*s)
dynamic_viscosity = 3.908e-7  # slugs/(ft * s)

# Calculate Reynolds number
reynolds_number = (airspeed * reference_chord_length * density) / dynamic_viscosity

print(f"Reynolds Number: {reynolds_number:.5e}")