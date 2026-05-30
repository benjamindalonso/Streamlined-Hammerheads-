import math
Ngear = 6 # Given in Raymers book page 370
Wlanding = 55000 # GTOW of our design
Vvertical = 20 # ft/s, based on carrier landing
g = 32.174 # ft/s^2
ShockEfficiency = 0.8 # Assumed shock absorption efficiency
Nmainshocks = 2
Nnoseshocks = 1

XNoseWheel = 7 # ft, distance from nose to nose wheel
XMainWheel = 32.8 # ft, distance from nose to main wheel
XAftCG = 30.765 # ft, distance from nose to aft CG
XForeCG = 28.514 # ft, distance from nose to fore CG
H = 8 # ft, distance between cg and ground
Na = XAftCG - XNoseWheel 
B = XMainWheel - XNoseWheel
Mf = XMainWheel - XForeCG
P = 259200 # Internal pressure in oleo shock in lbf/ft^2, based on Raymers book page 372

MSLmain = Wlanding * (Na / B) # Main gear static load
MSLnose = Wlanding * (Mf / B) # Nose gear static load
DBLnose = (10 * H * Wlanding) / (g * B) # Dynamic load on nose gear during landing

MainWheelOleoLoad = MSLmain / Nmainshocks # Load per main wheel shock
NoseWheelOleoLoad = (MSLnose + DBLnose) / Nnoseshocks # Load per nose wheel shock

MainOleoDiameter = 1.3 * math.sqrt((4 * MainWheelOleoLoad) / (math.pi * P))
NoseOleoDiameter = 1.3 * math.sqrt((4 * NoseWheelOleoLoad) / (math.pi * P))

# Calculate vertical shcock stroke required assuming the tire doesn't help with deflection (This eers on the safe side since the tire would help in real life)
S = (Vvertical ** 2) / (2 * g * ShockEfficiency * Ngear) # ft, shock stroke length
Smain = S + (1/12) # Add 1 additional inch as a safety factor (page 371 Raymer)
Snose = S + (2/12) # Add 2 additional inches becuase Raymer says nose gear is the same or slightly larger than main strut
print(f"Main Shock stroke length S = {Smain:.2f} ft")
print(f"Nose Shock stroke length S = {Snose:.2f} ft")

# Calculate length of oleo 
Lmain = Smain * 2.5 # Total oleo length is 2.5 times the stroke (Thisis the total length when there is no weight on wheels)
Lnose = Snose * 2.5
print(f"Main Oleo length L = {Lmain:.2f} ft")
print(f"Nose Oleo length L = {Lnose:.2f} ft")

print(f"Main Oleo diameter D = {MainOleoDiameter:.2f} ft")
print(f"Nose Oleo diameter D = {NoseOleoDiameter:.2f} ft")