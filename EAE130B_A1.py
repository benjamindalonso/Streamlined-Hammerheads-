g = 32.174
knot_to_fps = 1.68781

W = 46670          
x_cg = 28.3        
x_nose = 12.213     
x_main = 35.164    
H = 7.1           

n_main_wheels = 2
n_nose_wheels = 1
n_braked_wheels = 2

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