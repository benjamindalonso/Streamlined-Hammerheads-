import numpy as np
import matplotlib.pyplot as plt
import math

# Fixed parameters (Add more as needed to the bottom of this list)
MaxLD = 12
Range = 2000                  # nmi
Endurance = 0.333            # hr
Ct = 0.889                   # specific fuel concumption lbm/(lbf hr) 
V = 550                      # knots
e = 0.8 # Oswald efficiency factor (typical value for fighters)
HorizTailArea = 100
VertTailArea = 81
WetAreaFuse = 180
NumberOfEngines = 1  # Example number of engines
TurnRate = 0.1745 # in Rad/s based on rfp preference of 10 deg/s
g = 32.174 # ft/s^2 (gravitational acceleration)
Vturn = 750 # Velocity during the maneuver in feet per second
Cd0Turn = 0.00696 # Zero lift drag coefficient during the turn (I just used the openVSP value again)
AR = 2.5
K =  1/(math.pi*e*AR) # Induced drag factor
rhoTurn = 0.001267 # Slugs per ft^3 (Density at 20,000 ft - maneuvering altitude per rfp)
MidMissionFuelFraction =  0.906 # Fuel fraction half way through cruise portion of mission (This is based on the rfp requirements for maneuvering being done at "mid mission weight")
TakeoffFuelFraction = 0.99 # Fuel fraction from assignment2code
ClimbFuelFraction = 0.96 # Fuel fraction from assignment2code
ThrustReduction =  0.8 # Thrust reduction factor at cruise (due to altitude and speed)
WingArea = 610
E = 0.333         # hr loiter
Clmax = 1.5 # Maximum coefficient of lift (this will occur right at stall - max angle of attack)
ROC = 200 # Rate of climb ft/min
V_horizontal = 500 * 1.68781  # Climb airspeed knots to ft/s
V_horizontal_min = V_horizontal * 60  # Climb airspeed in ft/min
Climb_Gradient = ROC / V_horizontal_min # Climb gradient
Climb_Cd0 = 0.01696 # Climb drag coefficient (I just used the openVSP value again)
Ks = 1.8 # Stall speed factor 
rhoTropicalDay = 0.00219 # Air density at sea level on a tropical day in slugs/ft^3 (for launch constraint)
maxLandSpeed = 202.6 # Max landing speed in feet per second 
T_0 = 43000  # Example value for thrust per engine
TotalThrustInitialGuess = T_0 * NumberOfEngines
TOGW_Guess = 30000  # Initial guess for Takeoff Gross Weight in pounds
WingAreaGrid = list(range(200, 1400, 5))  # range of wing areas to analyze
Vend = 135 # Catipult end speed in knots with a 67,000 GTOW and a 210 CSV setting on the catipult 
Vwod = 0 # Wind speed over the deck in knots (Assumed 0 for worst case scenario)
Vthrust = 10 # Velocity added by engine thrust during catipult launch (Assumed to be 10 knots per Raymer page 136)
ClmaxTakeOff = 1.7 # Clmax at takeoff per slide 11 of preliminary sizing part 2

# Calculates weight of the crew
NumberOfPilots = 1 # Number of pilots
WeightOfPilot = 250 # Weight of pilot in lbf
WeightCrew = NumberOfPilots * WeightOfPilot
print("Crew Weight: " + str(WeightCrew) + " lbf ")

# Caluculates weight of the payload
MK83JDAM = 1000 # Weight of the MK-83 JDAM missle  in lbf
AIM9X = 190 # Weight of AIM-9X missle in lbf
MK83Quantity = 4 # Quantity of MK83 Missles
AIM9XQuantity = 2 # Quantity of AIM9x
Avionics = 2500 # Weight of avionics per RFP
WeightPayload = (MK83JDAM * MK83Quantity) + (AIM9X * AIM9XQuantity) + Avionics # Weight of paylod (which in this case is just armament)

# Engine weight calculation 
def Engine_Weight_Calculation(T_0):
    W_eng_dry = 0.521 * T_0**0.9
    W_eng_oil = 0.082 * T_0**0.65
    W_eng_rev = 0.034 * T_0
    W_eng_control = 0.26 * T_0**0.5
    W_eng_start = 9.33 * (W_eng_dry/1000) ** 1.078
    W_eng = W_eng_dry + W_eng_oil + W_eng_rev + W_eng_control + W_eng_start
    return W_eng

# Empty weight calculation
def Empty_Weight_Calculation(WingArea, HorizTailArea, VertTailArea, WetAreaFuse, TOGW, T_0 , NumberOfEngines):
    WingWeight = WingArea * 9
    HorizTailWeight = HorizTailArea * 4
    VertTailWeight = VertTailArea * 5.3
    FuseWeight = WetAreaFuse * 4.8
    LandingGearWeight = 0.045 * TOGW
    EWeight = Engine_Weight_Calculation(T_0)
    EngineWeight = EWeight * NumberOfEngines * 1.3
    AllElse = 0.17 * TOGW
    W_empty = WingWeight + HorizTailWeight + VertTailWeight + FuseWeight + LandingGearWeight + EngineWeight + AllElse
    return W_empty

# Weight Fraction Calculation 
def Weight_Fraction_Calculation(MaxLD, Range, Endurance, Ct, V):
    L_D = 0.94 * MaxLD
    W3_W2 = np.exp((-Range*Ct) / (V*L_D))  # cruise
    W4_W3 = np.exp((-Endurance*Ct) / (L_D))    # loiter/descent
    W1_W0 = 0.970   # engine start & takeoff
    W2_W1 = 0.96   # climb
    W5_W4 = 0.995   # landing
    W5_W0 = W5_W4 * W4_W3 * W3_W2 * W2_W1 * W1_W0
    Wf_W0 = (1 - W5_W0) * 1.06    # compute fuel fraction
    return Wf_W0



# This is where the constraints go (Basically what we did in assignment 3, don't put loops here)
def Climb_Constraint(Ks, K, Climb_Cd0, Clmax, Climb_Gradient):
    Climb_Intial = ((Ks**2*Climb_Cd0)/(Clmax))+((K*((Clmax)/Ks**2))+((Climb_Gradient)))
    Climb_Constraint = Climb_Intial*((1/.8)*(1/.99)) # Adjusting for fuel fractions and thrust reduction
    return Climb_Constraint

# Landing Constraint finds wing loading for landing weights lbf,speeds ft/s, density slug/ft^3
def Landing_Constraint(GTOW,landWeight,maxLandSpeed,CLmaxLand,density):
    S_land = 2*landWeight/(density*((maxLandSpeed/1.15)**2)*CLmaxLand)
    landWingLoading = GTOW/S_land
    return landWingLoading







#---------------------------------------------------------------------------
# Weight Inner Loop Calculation
def Weight_Inner_Loop(TOGW_Guess, WingArea, HorizTailArea, VertTailArea, WetAreaFuse, NumberOfEngines, WeightCrew, WeightPayload, T_0, err=1e-6, max_iter=200):
    W0_history = []
    delta = np.inf
    it = 0
    while delta > err and it < max_iter:
        # 1) Fuel Fraction
        Wf_W0 = Weight_Fraction_Calculation(MaxLD, Range, Endurance, Ct, V)
        FuelWeight = Wf_W0 * TOGW_Guess

        # 2) Empty Weight
        EmptyWeight = Empty_Weight_Calculation(WingArea, HorizTailArea, VertTailArea, WetAreaFuse, TOGW_Guess, T_0 , NumberOfEngines)

        # 3) New Gross Weight
        NewEmptyWeight = EmptyWeight + WeightCrew + WeightPayload + FuelWeight
        W0_history.append(NewEmptyWeight)

        # 4) convergence check
        delta = abs(NewEmptyWeight - TOGW_Guess) / max(abs(NewEmptyWeight), 1e-9)

        # 5) update
        TOGW_Guess = NewEmptyWeight
        it += 1

    converged = (delta <= err)
    return TOGW_Guess, converged, it, np.array(W0_history)
    
# Runs take off weight loop to come up with a weight estimation
W0, conv, it, hist = Weight_Inner_Loop(TOGW_Guess, WingArea,HorizTailArea, VertTailArea,WetAreaFuse, NumberOfEngines, WeightCrew, WeightPayload, T_0)

# Plots weight convergence
plt.figure()
plt.plot(hist, marker='o')
plt.xlabel("Iteration")
plt.ylabel("Takeoff Gross Weight W0")
plt.title("TOGW Convergence")
plt.grid(True)
plt.show()
#--------------------------------------------------------------------------------




# Here is where loops go

# Climb Loop
def outer_loop_climb_constraint(
    wing_area_grid=WingAreaGrid,           
    TOGW_guess_init=30000,
    T_total_guess_init=20000,
    tol_T_rel=1e-3,
    max_iter_T=60,
    relax=0.4):
   
    T_converged = []
    W0_converged = []
    
    for S_wing in wing_area_grid:
        T_total = T_total_guess_init
    
        
        for iter_outer in range(max_iter_T):
            T0_per_engine = T_total / NumberOfEngines
            
            W0, converged_inner, it_inner, _ = Weight_Inner_Loop(
                TOGW_Guess    = TOGW_guess_init,
                WingArea      = S_wing,
                HorizTailArea = HorizTailArea,
                VertTailArea  = VertTailArea,
                WetAreaFuse   = WetAreaFuse,
                NumberOfEngines = NumberOfEngines,
                WeightCrew    = WeightCrew,
                WeightPayload = WeightPayload,
                T_0           = T0_per_engine
            )
            
            if not converged_inner:
                print(" [inner loop did not converge]")
                break
            
            WS = W0 / S_wing
            
        
            climb_initial = (Ks**2 * Climb_Cd0 / Clmax) + \
                            (K * (Clmax / Ks**2)) + \
                            Climb_Gradient
            TW_req = climb_initial * (1 / 0.8) * (1 / 0.99)
            
            T_req = TW_req * W0
            
            rel_error = abs(T_req - T_total) / max(abs(T_total), 1e-6)
            if rel_error < tol_T_rel:
                T_converged.append(T_req)
                W0_converged.append(W0)
                break
            
            T_total = (1 - relax) * T_total + relax * T_req
        
        else:
            T_converged.append(T_req)
            W0_converged.append(W0)
    
    return np.array(T_converged), np.array(W0_converged)

# Runs climb constraint
T_climb, W0_climb = outer_loop_climb_constraint()


# Landing Loop
def outer_loop_landing_constraint(T_grid,TOGW_Guess_init,maxLandSpeed,CLmaxLand,density):
    tolerance = 10**(-6)
    S_converged = []
    for T in T_grid:
        S_wing = 1
        delta = 2*tolerance
        while delta > tolerance:
             W0, wconv, it_w, W0_hist = Weight_Inner_Loop(TOGW_Guess_init, S_wing, HorizTailArea, VertTailArea, WetAreaFuse, NumberOfEngines, WeightCrew, WeightPayload, T)
             fuelWeightLand = (W0 - (Empty_Weight_Calculation(S_wing,HorizTailArea,VertTailArea,WetAreaFuse,W0,T,NumberOfEngines) + WeightCrew + WeightPayload))*0.25
             landWeight =  fuelWeightLand + Empty_Weight_Calculation(S_wing,HorizTailArea,VertTailArea,WetAreaFuse,W0,T,NumberOfEngines) + WeightCrew + WeightPayload
             W_Sreq = Landing_Constraint(W0,landWeight,maxLandSpeed,CLmaxLand,density)
             Snew = W0/W_Sreq
             delta = abs(Snew - S_wing)
             S_wing = Snew
        S_converged.append(Snew)
    return np.array(S_converged)   

S_grid = np.linspace(300, 600, 7)
T_grid = np.linspace(0,70000,10)

# Runs landing constraint
S_converged_landing = outer_loop_landing_constraint(T_grid, TOGW_Guess, maxLandSpeed, Clmax, rhoTropicalDay)


# Launch loop
T_levels = np.linspace(0, 45000, 7)   # all the thrust values this loop will test
S_min_launch = []
for T in T_levels:
    T0 = T / NumberOfEngines
    W0, _, _, _ = Weight_Inner_Loop(
        TOGW_Guess    = TOGW_Guess,
        WingArea      = WingArea,
        HorizTailArea = HorizTailArea,
        VertTailArea  = VertTailArea,
        WetAreaFuse   = WetAreaFuse,
        NumberOfEngines = NumberOfEngines,
        WeightCrew    = WeightCrew,
        WeightPayload = WeightPayload,
        T_0           = T0
    )
    
    WS_max = 0.5 * rhoTropicalDay * ((Vend + Vwod + Vthrust)**2) * ClmaxTakeOff / 1.21
    S_min = W0 / WS_max
    S_min_launch.append(S_min)




# =============================================================================
#   FINAL CONSTRAINT DIAGRAM + SUMMARY OUTPUT
# ============================================================================

# ─── Print converged takeoff gross weight from initial estimation ───────────
print("\n" + "="*60)
print("INITIAL WEIGHT ESTIMATION (fixed wing area = {} ft²)".format(WingArea))
print("Converged Takeoff Gross Weight (W₀): {:.0f} lbf".format(W0))
print("  - Empty weight:       {:.0f} lbf".format(
    Empty_Weight_Calculation(WingArea, HorizTailArea, VertTailArea, WetAreaFuse, W0, T_0, NumberOfEngines)
))
print("  - Crew + Payload:     {:.0f} lbf".format(WeightCrew + WeightPayload))
print("  - Fuel weight:        {:.0f} lbf".format(W0 - (WeightCrew + WeightPayload + 
      Empty_Weight_Calculation(WingArea, HorizTailArea, VertTailArea, WetAreaFuse, W0, T_0, NumberOfEngines))))
print("  - Iterations to converge: {}".format(it))
print("="*60 + "\n")

# ─── Combined Constraint Diagram ────────────────────────────────────────────
plt.figure(figsize=(11, 7))

# Climb constraint: Thrust required vs Wing Area
plt.plot(WingAreaGrid, T_climb,
         marker='o', linestyle='-', linewidth=1.8, markersize=7,
         color='C0', label='Climb constraint (ROC + lapse adjustment)')

# Landing constraint: Thrust vs Required Wing Area (minimum S for given T)
plt.plot(S_converged_landing, T_grid,
         marker='s', linestyle='--', linewidth=1.8, markersize=7,
         color='C3', label='Landing constraint (max landing speed = {} ft/s)'.format(maxLandSpeed))

# Launch constraint: Thrust required vs wing area
plt.plot(S_min_launch, T_levels,
         marker='D', linestyle=':', linewidth=1.6, markersize=8,
         color='C4', label='Launch / Takeoff\n(min S to meet WS limit)')

# ADD NEW PLOTTING HERE FOLLOWING PREVIOUS FORMAT

# Formatting 
plt.xlabel("Wing Area S  (ft²)", fontsize=12)
plt.ylabel("Total Thrust Required  (lbf)", fontsize=12)
plt.title("Aircraft Sizing Constraint Diagram\n", fontsize=14, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=11, loc='upper right', framealpha=0.95, edgecolor='gray')
plt.tight_layout()
plt.show()

