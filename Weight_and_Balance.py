S_Design = 600 # ft^2
n_ultamate = 8 # Ultimate load factor
ARw = 2.028
W_Sweep = 38 # degrees
Taper_Ratio = .24
Root_Chord = 27.74 # feet
rho_30k = 0.000891 # slugs/ft^3
CruiseSpeed = 550 
Cruisefps = CruiseSpeed * 1.6878 
q_cruise = 0.5 * rho_30k * (Cruisefps**2)
W_design = 75000

S_ht = 77.27030 # ft^2
S_Vt = 86.21542 # ft^2
ARv = 2.88328
Bh = 19.36488 # ft
Fw = 10.6 # ft, fuselage width at wing intersection
Ht = .099 # Horizontal tail height above Fuselage reference
Hv = 15.192 # Vertical tail height above Fuselage reference
M = .8 # Mach number assumption
Tail_Sweep = 36.35714 # degrees
Lt = 44.027 # ft, tail moment arm
Sr = S_Vt*.3*.8 # ft^2, rudder area
L_fuse = 49.5 # ft, fuselage length
W_fuse = 11.41728 # ft, fuselage width
L_d = 8.16429 # ft, Length of intake
L_s = 3.2 # ft, length of straight intake duct guess
D_e = 43 # inches
D_e_ft = D_e / 12 # ft, equivalent diameter of engine exhaust




W_Fuselage = .499*(W_design**(.35))*(n_ultamate**(.25))*(L_fuse**(.5))*(W_fuse**(.849))*(W_design**(.685))
print(f"Fuselage Weight = {W_Fuselage:.2f} lbs")
W_Vt = .452*(((1+(Ht/Hv))**(.5))*((W_design*n_ultamate))**(.488)*(S_Vt**(.718)))*M**(.341)*Lt**(-1)*(1+Sr/S_Vt)**(.348)*ARv**(.223)*(1+Taper_Ratio)**(.25)*.805**(-.323) # lbs
print(f"Vertical Tail Weight = {W_Vt:.2f} lbs")
W_Ht = 3.316*(1+(Fw/Bh))**(-2)*((W_design*n_ultamate)/1000)**(.26)*(S_ht**(.806)) # lbs
print(f"Horizontal Tail Weight = {W_Ht:.2f} lbs")
W_wing = .036*(S_Design**(.758))*(ARw/(.6209))*(q_cruise**(.006))*(Taper_Ratio**(.04))*(((100*Root_Chord)/(.788))**(-.3))*((n_ultamate*W_design)**(.49)) # lbs
print(f"Wing Weight = {W_wing:.2f} lbs")
W_intake = 13.29*L_d**(.643)*(L_s/L_d)**(-.373)*D_e_ft 
print(f"Intake Weight = {W_intake:.2f} lbs")
W_Engine = 6422 # lbs