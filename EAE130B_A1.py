import numpy as np
import matplotlib.pyplot as plt
import math

g = 32.174
knot_to_fps = 1.68781

W = 54500          
x_cg =  29.659       
x_nose = 10    
x_main = 31.6590   
H = 7.787          

n_main_wheels = 2
n_nose_wheels = 2
n_braked_wheels = 2

brake_time = 25 # How many seconds brakes are applied (this guess needs to be refined/sourced properly)

V = 130 # Stall speed in knots             
a_brake = 10       

# geometry
B = x_main - x_nose           
Mf = x_cg - x_nose             
Ma = x_main - x_cg             

# static loads
main_load_total = W * Mf / B
nose_load_total = W * Ma / B

main_load_per_wheel = main_load_total / n_main_wheels
nose_load_per_wheel = nose_load_total / n_nose_wheels

# load on nose
dynamic_nose_load = (a_brake * H * W) / (g * B)
nose_load_braking = nose_load_total + dynamic_nose_load

# breaking kinetic energy
V_fps = V * knot_to_fps
KE_total = W * V_fps**2 / (2 * g)               
KE_per_wheel = KE_total / n_braked_wheels

KE_total_MJ = KE_total * 1.3558179483314 / 1e6
KE_per_wheel_MJ = KE_per_wheel * 1.3558179483314 / 1e6

# Convert to power units (energy per unit time)
avg_power_total = KE_total / brake_time                    # ft·lbf / s
avg_power_per_wheel = KE_per_wheel / brake_time            # ft·lbf / s  (per braked wheel)

# results
print("Wheelbase B =", round(B, 3), "ft")
print("Nose gear to CG Mf =", round(Mf, 3), "ft")
print("CG to main gear Ma =", round(Ma, 3), "ft")

print("\nstatic loads")
print("Total main gear load =", round(main_load_total, 1), "lbf")
print("Total nose gear load =", round(nose_load_total, 1), "lbf")
print("Main wheel load =", round(main_load_per_wheel, 1), "lbf per wheel")
print("Nose wheel load =", round(nose_load_per_wheel, 1), "lbf per wheel")

print("\nbraking load")
print("Additional nose load during braking =", round(dynamic_nose_load, 1), "lbf")
print("Total nose load during braking =", round(nose_load_braking, 1), "lbf")

print("\nbraking kinetic energy")
print("Speed =", round(V_fps, 2), "ft/s")
print("Total KE =", round(KE_total, 1), "ft*lbf")
print("KE per braked wheel =", round(KE_per_wheel, 1), "ft*lbf")
print("Total KE =", round(KE_total_MJ, 3), "MJ")
print("KE per braked wheel =", round(KE_per_wheel_MJ, 3), "MJ")
print("Average power per wheel =", round(avg_power_per_wheel, 1), "ft*lbf/s")

# Tire Sizing

main_wheel_diameter = 1.59 * (main_load_per_wheel ** 0.302) # Fighter/trainer relation from raymer page 358
main_wheel_width = 0.0980 * (main_load_per_wheel ** 0.467) # Fighter/trainer relation from raymer page 358
nose_wheel_diameter = 0.7 * main_wheel_diameter # Nose wheel diameter is %60-%100 of main wheel diameter according to slide 47 of lecture 1
nose_wheel_width = 0.7 * main_wheel_width # Lecture doesn't specify how to choose nose tire width, here we assume it's similar to the diameter

print("Main wheel diameter =", round(main_wheel_diameter, 1), "in")
print("Main wheel width =", round(main_wheel_width, 1), "in")
print("Nose wheel diameter =", round(nose_wheel_diameter, 1), "in")
print("Nose wheel width =", round(nose_wheel_width, 1), "in")

# Tire sizes need to be compared to graph on page 50 of lecture 1 to determine if the wheels are big enough for the required brakes (see average power per wheel)











# Landing gear geometry checks 
x_cg_fully_loaded = 29.659
z_cg_fully_loaded = -0.377

x_cg_nofuel_noarmamants = 27.347
z_cg_nofuel_noarmamants = -0.289

x_furthest_aft_and_lowest_point = 47.786
z_furthest_aft_and_lowest_point = -3.557

x_main_gear = 31.6590 # Feet from nose of the center of the main gear wheels
z_main_gear = -8.0760 # Feet below the nose of the bottom of the main gear wheels

# Longitudinal Tip Over and Ground Clearance (these two tests basically determine where the main gear has to be since it's the intersection of two lines)
x_start1 = x_cg_fully_loaded
z_start1 = z_cg_fully_loaded
x_start2 = x_furthest_aft_and_lowest_point
z_start2 = z_furthest_aft_and_lowest_point

# Angles
theta1 = np.radians(-75)   # from fully loaded CG
theta2 = np.radians(16)    # from rear-most point

# Direction vectors
dx1 = np.cos(theta1)
dz1 = np.sin(theta1)
dx2 = np.cos(theta2)
dz2 = np.sin(theta2)

# Create the plot
plt.figure(figsize=(14, 9))

# Plot the two starting points
plt.plot(x_start1, z_start1, 'bo', markersize=10, label='Fully Loaded CG')
plt.plot(x_start2, z_start2, 'ro', markersize=10, label='Rear-most Point')

# Plot Line 1: from CG at -75° (extended in both directions)
t_range = np.linspace(-40, 150, 600)
x_line1 = x_start1 + t_range * dx1
z_line1 = z_start1 + t_range * dz1
plt.plot(x_line1, z_line1, 'b-', linewidth=2.8, label='Line from CG at -75°')

# Plot Line 2: from rear at +16°
s_range = np.linspace(-150, 100, 600)
x_line2 = x_start2 + s_range * dx2
z_line2 = z_start2 + s_range * dz2
plt.plot(x_line2, z_line2, 'r-', linewidth=2.8, label='Line from Rear at +16°')

plt.xlabel('X (ft) — Forward')
plt.ylabel('Z (ft) — Vertical (negative = down)')
plt.title('Graphical Intersection Selection\n'
          'Click on the plot where the two lines cross')

plt.grid(True, alpha=0.4)
plt.axis('equal')
plt.legend()

# Clicking anywhere on the plot prints the coordinates
def onclick(event):
    if event.xdata is not None and event.ydata is not None:
        print(f"Clicked →  x = {event.xdata:.4f} ft,   z = {event.ydata:.4f} ft")

fig = plt.gcf()
cid = fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()


# Lateral Tip Over Test

# adjust these two values to acheive desired tip over angle (these do not effect ground clearnace or longitudinal tip over)
y_main_gear = 7 # Feet from centerline to middle of main gear wheel
x_nose_gear = 7 # Feet from nose to center of nose gear


crit_angle = 55 # Degrees
crit_angle_rad = np.radians(crit_angle)
h = abs(z_main_gear - z_cg_nofuel_noarmamants) # height of cg above bottom point on main wheels
nose_main_slope = math.pi/2 - math.atan((x_main_gear - x_nose_gear) / y_main_gear) # slope of line from nose gear to mian gear when viewed from above
line_seperation = (x_cg_nofuel_noarmamants - x_nose_gear) * math.sin(nose_main_slope) # separation between line from nose gear to main gear and parallel line through the CG
lateral_tip_over_angle = math.atan(h / line_seperation)
lateral_tip_over_angle_deg = np.degrees(lateral_tip_over_angle)
print("\nLateral Tip Over Test")
print("Lateral tip over angle =", round(lateral_tip_over_angle_deg, 2), "degrees")
if lateral_tip_over_angle < crit_angle_rad:
    print("Lateral tip over test PASSED")
else:
    print("Lateral tip over test FAILED")


