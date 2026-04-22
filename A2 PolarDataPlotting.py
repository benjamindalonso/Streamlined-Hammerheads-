# NACA 65-410 XFOIL polar data (Mach=0.0, Re=25.32e6, Ncrit=9.0)
# alpha (degrees) and sectional lift coefficient (unitless)
import os
import matplotlib.pyplot as plt

# ====================== NACA 65-410 Landing AIRFOIL PLOT ======================

x = [
    0.932823, 0.900766, 0.868495, 0.836040, 0.803281, 0.789288,
    0.744737, 0.694687, 0.644567, 0.594367, 0.544127, 0.493837,
    0.443517, 0.393197, 0.342867, 0.292557, 0.242267, 0.192007,
    0.141817, 0.087697, 0.061381, 0.039089, 0.017382, 0.006968,
    0.003117, 0.001365, 0.000000, 0.007777, 0.010647, 0.016072,
    0.028989, 0.054233, 0.079118, 0.090009, 0.095450, 0.112243,
    0.145857, 0.195667, 0.245407, 0.295117, 0.344807, 0.394477,
    0.444157, 0.493837, 0.543547, 0.593307, 0.643107, 0.692987,
    0.742937, 0.762762, 0.772677, 0.787532, 0.817096, 0.857234,
    0.896400, 0.932823
]

y = [
    -0.108781, -0.069650, -0.030758, 0.008054, 0.046612, 0.060878,
    0.069239, 0.076509, 0.082929, 0.088399, 0.092719, 0.095699,
    0.097049, 0.096899, 0.095349, 0.092539, 0.088419, 0.082829,
    0.075579, 0.065178, 0.052547, 0.039461, 0.024690, 0.016071,
    0.011820, 0.009279, 0.000000, -0.004759, -0.005233, -0.005466,
    -0.004563, -0.001370, 0.002658, 0.004580, 0.005568, 0.004375,
    0.002379, 0.000049, -0.001581, -0.002621, -0.003111, -0.003021,
    -0.002211, -0.000541, 0.002119, 0.005479, 0.009309, 0.013409,
    0.017599, 0.019215, 0.020007, 0.009351, -0.012158, -0.042404,
    -0.073950, -0.108781
]

# Plot the airfoil
plt.figure(figsize=(12, 5))
plt.plot(x, y, 'b-', linewidth=2.5, label='Flaps -40 Deg, Slats -15 Deg')
plt.fill(x, y, color='lightblue', alpha=0.4)

plt.xlabel('x/c  (chord fraction)', fontsize=12)
plt.ylabel('y/c', fontsize=12)
plt.title('NACA 65-410 Airfoil Geometry\nLanding Configuration', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.axis('equal')
plt.legend(fontsize=11)
plt.tight_layout()
plt.show()

# ====================== NACA 65-410 Takeoff AIRFOIL PLOT ======================

x = [
    0.978878, 0.935371, 0.891745, 0.847972, 0.804001, 0.789591,
    0.748445, 0.698395, 0.648275, 0.598075, 0.547835, 0.497545,
    0.447225, 0.396905, 0.346575, 0.296265, 0.245975, 0.195715,
    0.145525, 0.092269, 0.067328, 0.043606, 0.020261, 0.008875,
    0.004534, 0.002486, 0.000000, 0.007139, 0.009930, 0.015286,
    0.028217, 0.053661, 0.078852, 0.091144, 0.097288, 0.114704,
    0.149565, 0.199375, 0.249115, 0.298825, 0.348515, 0.398185,
    0.447865, 0.497545, 0.547255, 0.597015, 0.646815, 0.696695,
    0.746645, 0.769168, 0.780432, 0.799341, 0.837084, 0.885146,
    0.932740, 0.978878
]

y = [
    -0.058715, -0.032909, -0.007399, 0.017972, 0.043000, 0.050704,
    0.057534, 0.064804, 0.071224, 0.076694, 0.081014, 0.083994,
    0.085344, 0.085194, 0.083644, 0.080834, 0.076714, 0.071124,
    0.063874, 0.053756, 0.044675, 0.034404, 0.022388, 0.015102,
    0.011352, 0.009044, 0.000000, -0.005672, -0.006491, -0.007384,
    -0.008061, -0.007969, -0.007004, -0.006353, -0.005995, -0.007248,
    -0.009326, -0.011656, -0.013286, -0.014326, -0.014816, -0.014726,
    -0.013916, -0.012246, -0.009586, -0.006226, -0.002396, 0.001704,
    0.005894, 0.007726, 0.008621, 0.003334, -0.007500, -0.022195,
    -0.038442, -0.058715
]

# Plot the airfoil
plt.figure(figsize=(12, 5))
plt.plot(x, y, 'r-', linewidth=2.5, label='Flaps -20 Deg, Slats -8 Deg')
plt.fill(x, y, color='lightcoral', alpha=0.4)

plt.xlabel('x/c  (chord fraction)', fontsize=12)
plt.ylabel('y/c', fontsize=12)
plt.title('NACA 65-410 Airfoil Geometry\nTakeoff Configuration', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.axis('equal')
plt.legend(fontsize=11)
plt.tight_layout()
plt.show()


# ====================== NACA 65-410 Clean (No high lift devices) AIRFOIL PLOT ======================

x = [
    1.00000, 0.95029, 0.90057, 0.85076, 0.80088, 0.75090,
    0.70085, 0.65073, 0.60053, 0.55029, 0.50000, 0.44968,
    0.39936, 0.34903, 0.29872, 0.24843, 0.19817, 0.14798,
    0.09788, 0.07289, 0.04797, 0.02318, 0.01089, 0.00607,
    0.00372, 0.00000, 0.00628, 0.00893, 0.01411, 0.02682,
    0.05203, 0.07711, 0.10212, 0.15202, 0.20183, 0.25157,
    0.30128, 0.35097, 0.40064, 0.45032, 0.50000, 0.54971,
    0.59947, 0.64927, 0.69915, 0.74910, 0.79912, 0.84924,
    0.89943, 0.94971, 1.00000
]

y = [
    0.00000, 0.00937, 0.01842, 0.02729, 0.03577, 0.04372,
    0.05099, 0.05741, 0.06288, 0.06720, 0.07018, 0.07153,
    0.07138, 0.06983, 0.06702, 0.06290, 0.05731, 0.05006,
    0.04067, 0.03487, 0.02800, 0.01935, 0.01372, 0.01061,
    0.00861, 0.00000, -0.00661, -0.00781, -0.00944, -0.01191,
    -0.01536, -0.01791, -0.01999, -0.02314, -0.02547, -0.02710,
    -0.02814, -0.02863, -0.02854, -0.02773, -0.02606, -0.02340,
    -0.02004, -0.01621, -0.01211, -0.00792, -0.00393, -0.00037,
    0.00226, 0.00327, 0.00000
]

# Plot the airfoil
plt.figure(figsize=(12, 5))
plt.plot(x, y, 'g-', linewidth=2.5, label='Flaps 0 Deg, Slats 0 Deg')
plt.fill(x, y, color='lightgreen', alpha=0.4)

plt.xlabel('x/c  (chord fraction)', fontsize=12)
plt.ylabel('y/c', fontsize=12)
plt.title('NACA 65-410 Airfoil Geometry\nClean Configuration', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.axis('equal')
plt.legend(fontsize=11)
plt.tight_layout()
plt.show()



#===========================Polar Data Plotting============================
alphaLanding = [
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
    5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
    10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5,
    15.0, 15.5, 16.0
]

clLanding = [
    1.6319, 1.6751, 1.7178, 1.7598, 1.8014, 1.8425, 1.8829, 1.9227, 1.9617, 1.9997,
    2.0369, 2.0732, 2.1087, 2.1432, 2.1762, 2.2074, 2.2377, 2.2621, 2.2918, 2.3147,
    2.3367, 2.3578, 2.3709, 2.3848, 2.3944, 2.4051, 2.4120, 2.4154, 2.4163, 2.4145,
    2.4081, 2.4021, 2.3956 
]

alphaTakeoff = [
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
    5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
    10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5,
    15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5,
    20.0
]

clTakeoff = [
    1.2409, 1.2858, 1.3318, 1.3772, 1.4222, 1.4671, 1.5121, 1.5563, 1.6006, 1.6444,
    1.6876, 1.7299, 1.7719, 1.8199, 1.8603, 1.8993, 1.9376, 1.9755, 2.0113, 2.0460,
    2.0796, 2.1092, 2.1380, 2.1674, 2.1931, 2.2148, 2.2367, 2.2548, 2.2739, 2.2855,
    2.2996, 2.3043, 2.3088, 2.3097, 2.3065, 2.3016, 2.2945, 2.2833, 2.2680, 2.2496,
    2.2260
]


alphaClean = [
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
    5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.5,
    10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5,
    14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5,
    18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5,
    22.0, 22.5, 23.0
]

clClean = [
    0.1140, 0.1746, 0.2361, 0.2983, 0.3618, 0.4248, 0.4868, 0.5470, 0.6062, 0.6651,
    0.7826, 0.8410, 0.8991, 0.9569, 1.0142, 1.0712, 1.1275, 1.2363,
    1.2915, 1.3458, 1.3991, 1.4517, 1.5027, 1.5522, 1.6000, 1.6463,
    1.6891, 1.7299, 1.7677, 1.7969, 1.8140, 1.8306, 1.8450, 1.8595,
    1.8709, 1.8775, 1.8793, 1.8754, 1.8669, 1.8519, 1.8292, 1.7958,
    1.7521, 1.6905, 1.6012
]


# Find max CL and corresponding alpha for each configuration
max_cl_landing = max(clLanding)
max_alpha_landing = alphaLanding[clLanding.index(max_cl_landing)]

max_cl_takeoff = max(clTakeoff)
max_alpha_takeoff = alphaTakeoff[clTakeoff.index(max_cl_takeoff)]

max_cl_clean = max(clClean)
max_alpha_clean = alphaClean[clClean.index(max_cl_clean)]

# Create the plot
plt.figure(figsize=(14, 11))

plt.plot(alphaLanding, clLanding, 'b-', linewidth=2.5, label='Landing MAC (Flaps -40 Deg, Slats -15 Deg)')
plt.plot(alphaTakeoff, clTakeoff, 'r-', linewidth=2.5, label='Takeoff MAC (Flaps -20 Deg, Slats -8 Deg)')
plt.plot(alphaClean, clClean, 'g-', linewidth=2.5, label='Clean Tip (Flaps 0 Deg, Slats 0 Deg)')

# Mark the maximum CL points
plt.plot(max_alpha_landing, max_cl_landing, 'bo', markersize=8)
plt.plot(max_alpha_takeoff, max_cl_takeoff, 'ro', markersize=8)
plt.plot(max_alpha_clean, max_cl_clean, 'go', markersize=8)

# Annotate the max CL values
plt.annotate(f'Max Cl = {max_cl_landing:.2f}',
             xy=(max_alpha_landing, max_cl_landing),
             xytext=(10, 15), textcoords='offset points',
             arrowprops=dict(arrowstyle='->', color='blue'),
             fontsize=14, color='blue')

plt.annotate(f'Max Cl = {max_cl_takeoff:.2f}',
             xy=(max_alpha_takeoff, max_cl_takeoff),
             xytext=(10, 15), textcoords='offset points',
             arrowprops=dict(arrowstyle='->', color='red'),
             fontsize=14, color='red')

plt.annotate(f'Max Cl = {max_cl_clean:.2f}',
             xy=(max_alpha_clean, max_cl_clean),
                xytext=(10, 15), textcoords='offset points',
                arrowprops=dict(arrowstyle='->', color='green'),
                fontsize=14, color='green')


plt.xlabel('Angle of Attack α (degrees)', fontsize=14)
plt.ylabel('Sectional Lift Coefficient Cl', fontsize=14)
plt.title('NACA 65-410 Lift Curves\nVarious Configurations', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=14)
plt.tight_layout()
plt.show()   

# Print max values
print(f"Landing Config  →  Max CL = {max_cl_landing:.4f}  at  α = {max_alpha_landing}°")
print(f"Takeoff Config  →  Max CL = {max_cl_takeoff:.4f}  at  α = {max_alpha_takeoff}°")
print(f"Clean Config    →  Max CL = {max_cl_clean:.4f}  at  α = {max_alpha_clean}°")