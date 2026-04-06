g = 32.174
knot_to_fps = 1.68781

W = 46670          
x_cg = 28.3        
x_nose = 12.213     
x_main = 35.164    
H = 7.1           

n_main_wheels = 2
n_nose_wheels = 2
n_braked_wheels = 2

brake_time = 25 # How many seconds brakes are applied (this guess needs to be refined/sourced properly)

V = 140             
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