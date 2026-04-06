import numpy as np
import matplotlib.pyplot as plt
import math

# Parameters (from prior assignments and pdr)
W_start = 58619.0      # Weight at start of cruise from pdr (lbs)
S_ref = 600.0          # Wing area (ft^2)
CD0 = 0.00766           # Zero-lift drag coefficient
AR = 2.028               # Aspect Ratio
e = 0.8             # Oswald efficiency factor
k = 1.0 / (np.pi * AR * e) # Induced drag factor

V_inf_knots = 550    # Cruise speed in knots
V_inf = V_inf_knots * 1.68781 # Convert knots to ft/s
rho = 0.000891          # Air density at cruise altitude (slugs/ft^3)
c_t_hr = 0.8           # Thrust specific fuel consumption (1/hr)
c_t = c_t_hr / 3600.0  # Convert TSFC to 1/seconds

Total_Cruise_Range_nmi = 1400.0 # Total range minus climb distance
Total_Cruise_Range_ft = Total_Cruise_Range_nmi * 6076.12 # Convert to ft

# Segments to test for the error plot
segment_counts = [1, 4, 10, 40] 

# Exponential Cruise Discretization
def calculate_cruise_profile(segments):
    dR = Total_Cruise_Range_ft / segments
    
    weights = [W_start]
    ranges = [0.0]
    weight_fractions = [1.0]
    
    W_current = W_start
    
    for i in range(segments):
        # Calculate CL for weight
        CL_i = (2.0 * W_current) / (rho * (V_inf**2) * S_ref)
        
        # Calculate L/D for CL
        LD_i = CL_i / (CD0 + k * (CL_i**2))
        
        # Calculate weight at end of this segment (Exponential Approach)
        exponent = - (dR * c_t) / (V_inf * LD_i)
        W_next = W_current * np.exp(exponent)
        
        # Store
        weights.append(W_next)
        ranges.append((i + 1) * dR / 6076.12) # Store range in nmi
        weight_fractions.append(W_next / W_start)
        
        W_current = W_next
        
    return ranges, weight_fractions, weights

# Plots and Data
plt.figure(figsize=(12, 5))

# Weight Fraction vs Range plot
plt.subplot(1, 2, 1)
baseline_ranges, baseline_fractions, _ = calculate_cruise_profile(1000) # "Analytical" baseline (high segment count)
plt.plot(baseline_ranges, baseline_fractions, 'k--', label='Analytical (1000 seg)')

colors = ['g', 'r', 'c', 'm']
results = {}

for idx, segs in enumerate(segment_counts):
    rngs, fracs, wts = calculate_cruise_profile(segs)
    results[segs] = (rngs, fracs, wts)
    plt.plot(rngs, fracs, marker='o', color=colors[idx], label=f'{segs} segments')

plt.title('Weight Fraction vs Range')
plt.xlabel('Range [nmi]')
plt.ylabel('Weight Fraction ($W / W_{start}$)')
plt.legend()
plt.grid(True)

# Weight Error vs Range plot
plt.subplot(1, 2, 2)
# Calculate error relative to baseline
for idx, segs in enumerate(segment_counts):
    rngs, fracs, _ = results[segs]
    
    # Interpolate baseline for error
    baseline_interp = np.interp(rngs, baseline_ranges, baseline_fractions)
    error_percentage = np.abs(fracs - baseline_interp) / baseline_interp * 100.0
    
    plt.plot(rngs, error_percentage, marker='o', color=colors[idx], label=f'{segs} segments')

plt.title('Weight Error vs Range')
plt.xlabel('Range [nmi]')
plt.ylabel('Weight Error [%]')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# Print final values for tables
print(f"--- Cruise Segment Results ---")
final_ranges, final_fractions, final_weights = results[40] # Using 40 segments for high accuracy
print(f"Initial Cruise Weight: {final_weights[0]:.2f} lbs")
print(f"Final Cruise Weight:   {final_weights[-1]:.2f} lbs")
print(f"Fuel Burned in Cruise: {final_weights[0] - final_weights[-1]:.2f} lbs")
print(f"Cruise Weight Fraction (W_end / W_start): {final_fractions[-1]:.4f}")

# Optm cruise climb density profile
CL_opt = np.sqrt(CD0 / k)
LD_max = CL_opt / (CD0 + k * CL_opt**2)

print(f"\n--- Optimal Cruise Parameters ---")
print(f"Optimal CL for (L/D)max: {CL_opt:.4f}")
print(f"Maximum L/D: {LD_max:.2f}")

final_ranges, final_fractions, final_weights = results[40]
weights_array = np.array(final_weights)

# From W = 0.5 * rho * V^2 * S * CL_opt
required_densities = (2.0 * weights_array) / ((V_inf**2) * S_ref * CL_opt)

# Required Density vs Range Plot
plt.figure(figsize=(8, 5))
plt.plot(final_ranges, required_densities, 'b-', linewidth=2, label='Required Density $\\rho$')

# Formatting the plot
plt.title('Required Density Profile for Optimal Cruise Climb at $(L/D)_{max}$')
plt.xlabel('Range [nmi]')
plt.ylabel('Air Density $\\rho$ [slugs/ft$^3$]')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.show()

# Alt profile plot (cruise climb)

# Standard Atmo Function
def density_to_altitude(rho):
    """slugs/ft^3 to feet based on Std Atm"""
    rho_sl = 0.0023769        # Density at sea level
    rho_strat = 0.00070612    # Density at tropopause (36,089.2 ft)
    h_strat = 36089.2         # Altitude of tropopause
    
    if rho > rho_strat:
        # Troposphere logic (Standard Temperature Lapse Rate)
        h = (1 - (rho / rho_sl)**(1 / 4.25588)) / (6.8755856e-6)
    else:
        # Stratosphere logic (Isothermal Layer)
        h = h_strat - 20806.2 * np.log(rho / rho_strat)
    return h

# Convert required_densities array to altitudes
altitudes = [density_to_altitude(rho) for rho in required_densities]