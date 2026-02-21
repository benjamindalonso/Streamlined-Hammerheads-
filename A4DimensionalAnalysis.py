import numpy as np
import matplotlib.pyplot as plt
import math

# Our Chosen Design Point
S_Design = 665 # Ft^2
T_Design = 43000 # lbf

# Fighter T vs S comparison

S_F22 = 840 # ft^2
T_W_F22 = 1.4 # F-22 Raptor has a thrust-to-weight ratio of about 1.4 at takeoff
TOGW_F22 = 60000 # lbf
T_F22 = T_W_F22 * TOGW_F22 # lbf

S_F35 = 460.05
T_F35 = 43000

S_DRM = 45.7*10.764 # Wing Area for Dassalt Rafale M (m^2 ---> ft^2)
T_DRM = 16860*2 # Thrust for Dassault Rafale M (lbf)

S_F18 = 500 # F-18 wing area(ft^2)
T_F18 = 44000 # F-18 thrust (lbf)

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
TurnRate = 0.139626 # in Rad/s based on rfp minimum of 8 deg/s
TurnRate_10deg = 0.174533 # in Rad/s based on rfp minimum of 10 deg/s
g = 32.174 # ft/s^2 (gravitational acceleration)
Vturn = 750 # Velocity during the maneuver in feet per second
Cd0Turn = 0.00696 # Zero lift drag coefficient during the turn (I just used the openVSP value again)
AR = 2.5
K =  1/(math.pi*e*AR) # Induced drag factor
rhoTurn = 0.001267 # Slugs per ft^3 (Density at 20,000 ft - maneuvering altitude per rfp)
MidMissionFuelFraction =  0.83 # Fuel fraction halfway through cruise portion of mission (This is based on the rfp requirements for maneuvering being done at "mid mission weight")
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
TOGW_Guess = 50000  # Initial guess for Takeoff Gross Weight in pounds
WingAreaGrid = list(range(100, 1200, 5))  # range of wing areas to analyze
WingAreaGridManeuver = list(range(475,1200,5)) # The thrust required blows up to infinity at low wing areas, this fixes that
Vend = 135 # Catapult end speed in knots with a 67,000 GTOW and a 210 CSV setting on the catapult 
Vwod = 0 # Wind speed over the deck in knots (Assumed 0 for worst case scenario)
Vthrust = 10 # Velocity added by engine thrust during catapult launch (Assumed to be 10 knots per Raymer page 136)
Vendfps = 1.6878 * Vend # This converts knots to feet per second
Vwodfps = 1.6878 * Vwod # This converts knots to feet per second
Vthrustfps = 1.6878 * Vthrust # This converts knots to feet per second
ClmaxTakeOff = 1.7 # Clmax at takeoff per slide 11 of preliminary sizing part 2
rho_30k = 0.000889 # slug/ft^3 
V_fts = V * 1.68781 # ft/s 
q_cruise = 0.5 * rho_30k * V_fts**2
CD0_cruise = 0.00696 
a_dash = 994.0 # ft/s
M_dash = 2.0
V_dash = M_dash * a_dash
q_dash = 0.5 * rho_30k * V_dash**2
CD0_dash = CD0_cruise

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
    Climb_Constraint = Climb_Intial*((1/ThrustReduction)*(1/TakeoffFuelFraction)) # Adjusting for fuel fractions and thrust reduction
    return Climb_Constraint

# Landing Constraint finds wing loading for landing weights lbf,speeds ft/s, density slug/ft^3
def Landing_Constraint(GTOW,landWeight,maxLandSpeed,CLmaxLand,density):
    S_land = 2*landWeight/(density*((maxLandSpeed/1.15)**2)*CLmaxLand)
    landWingLoading = GTOW/S_land
    return landWingLoading


# Maneuvering Constraint (from A3)
def Sustained_Turn_Constraint(TurnRate, g, Vturn, Cd0, K, Wing_Loading, rho, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction):
    n = math.sqrt(((TurnRate * Vturn) / g)**2 + 1)
    q = 0.5 * rho * (Vturn**2)
    Wing_LoadingTurn = Wing_Loading * TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction 
    
    # Required T/W at maneuver alt
    TW_turn = (q * Cd0 / Wing_LoadingTurn) + (K * (n**2) * Wing_LoadingTurn / q)
    
    # Correct back to Sea Level Static T/W (T0/W0)
    TW_SLS = TW_turn * (TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction / ThrustReduction)
    return TW_SLS


#---------------------------------------------------------------------------
# Weight Inner Loop Calculation
def Weight_Inner_Loop(TOGW_Guess, WingArea, HorizTailArea, VertTailArea, WetAreaFuse, NumberOfEngines, WeightCrew, WeightPayload, T_0, err=1e-6, max_iter=1000):
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
    maxLandSpeed = maxLandSpeed + 25.3 # Adds 15 kts wind over deck in ft/s
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
S_converged_landing = outer_loop_landing_constraint(T_grid, TOGW_Guess, maxLandSpeed, 2, rhoTropicalDay)

# Launch loop
T_levels_launch = np.linspace(20000, 75000, 10)   # adjust range as needed
S_min_launch = []
T_for_plot = []  # to plot T vs min S

for T_total in T_levels_launch:
    T0 = T_total / NumberOfEngines
    WS_max = 0.5 * rhoTropicalDay * ((Vendfps + Vwodfps + Vthrustfps)**2) * ClmaxTakeOff / 1.21

    
    min_S_found = np.inf
    feasible = False
    
    for S_wing in WingAreaGrid:
        W0, converged, _, _ = Weight_Inner_Loop(
            TOGW_Guess    = TOGW_Guess,
            WingArea      = S_wing,          # ← use grid values here
            HorizTailArea = HorizTailArea,
            VertTailArea  = VertTailArea,
            WetAreaFuse   = WetAreaFuse,
            NumberOfEngines = NumberOfEngines,
            WeightCrew    = WeightCrew,
            WeightPayload = WeightPayload,
            T_0           = T0,
            err=1e-5,
            max_iter=1000
        )
        
    
        
        WS = W0 / S_wing
        
        
        if WS <= WS_max:
            feasible = True
            if S_wing < min_S_found:
                min_S_found = S_wing
    
    if feasible:
        S_min_launch.append(min_S_found)
        T_for_plot.append(T_total)
    else:
        S_min_launch.append(np.nan)
        T_for_plot.append(T_total)

# Maneuver loop
def outer_loop_maneuver_constraint(
        WingAreaGrid,
        TOGW_Guess,
        TotalThrustInitialGuess,
        TurnRate_input,
        tol_T_rel=1e-2,
        max_iter_T=100,
        relax=1):

    T_total_converged = []
    W0_converged = []

    for WingArea in WingAreaGridManeuver:

        T_total = TotalThrustInitialGuess

        for k in range(max_iter_T):

            # Thrust per engine
            T_0 = T_total / NumberOfEngines
        

            # Inner weight loop
            Final_TOGW, Converged, Iterations, _ = Weight_Inner_Loop(
                TOGW_Guess,
                WingArea,
                HorizTailArea,
                VertTailArea,
                WetAreaFuse,
                NumberOfEngines,
                WeightCrew,
                WeightPayload,
                T_0
            )
            WingLoading = Final_TOGW / WingArea
            
        
           
            # Maneuver T/W requirement
            TW_maneuver = Sustained_Turn_Constraint(
                TurnRate_input,
                g,
                Vturn,
                Cd0Turn,
                K,
                WingLoading,
                rhoTurn,
                MidMissionFuelFraction,
                TakeoffFuelFraction,
                ClimbFuelFraction,
                ThrustReduction
            )
            T_required = TW_maneuver * Final_TOGW
            # Checking the convergence
            rel_error = abs(T_required - T_total) / max(abs(T_total), 1e-9)

            if rel_error < tol_T_rel:
                T_total = T_required
                break
            # Relaxation update
            T_total = (1 - relax) * T_total + relax * T_required

        T_total_converged.append(T_total)
        W0_converged.append(Final_TOGW)

    return np.array(T_total_converged), np.array(W0_converged)

T_maneuver, W0_curve = outer_loop_maneuver_constraint(
    WingAreaGrid,
    TOGW_Guess,
    TotalThrustInitialGuess, TurnRate
)

T_maneuver_10, _ = outer_loop_maneuver_constraint(
    WingAreaGrid,
    TOGW_Guess,
    TotalThrustInitialGuess,
    TurnRate_10deg
)

def cruise_TW_req_SLS(WingLoading_W0S, q_cruise, CD0_cruise,
                      FuelFrac_cruise, ThrustReduction_cruise):
    WS_cruise = WingLoading_W0S * FuelFrac_cruise
    TW_alt = (q_cruise * CD0_cruise) / WS_cruise + (K * WS_cruise) / q_cruise
    TW_SLS = TW_alt * (FuelFrac_cruise / ThrustReduction_cruise)
    return TW_SLS

def outer_loop_cruise_constraint(
    wing_area_grid=WingAreaGrid,
    TOGW_guess_init=30000,
    T_total_guess_init=20000,
    CD0_cruise=0.00696,    
    q_cruise=300.0,
    FuelFrac_cruise=0.90, 
    ThrustReduction_cruise=0.8,
    tol_T_rel=1e-3,
    max_iter_T=60,
    relax=0.4
):
    T_converged = []
    W0_converged = []

    for S_wing in wing_area_grid:
        T_total = T_total_guess_init

        for iter_outer in range(max_iter_T):
            T0_per_engine = T_total / NumberOfEngines

            W0, converged_inner, it_inner, _ = Weight_Inner_Loop(
                TOGW_Guess       = TOGW_guess_init,
                WingArea         = S_wing,
                HorizTailArea    = HorizTailArea,
                VertTailArea     = VertTailArea,
                WetAreaFuse      = WetAreaFuse,
                NumberOfEngines  = NumberOfEngines,
                WeightCrew       = WeightCrew,
                WeightPayload    = WeightPayload,
                T_0              = T0_per_engine
            )

            if not converged_inner:
                print(f"[inner loop did not converge at S={S_wing}]")
                break

            WS_W0 = W0 / S_wing  # W0/S

            TW_req_SLS = cruise_TW_req_SLS(
                WingLoading_W0S=WS_W0,
                q_cruise=q_cruise,
                CD0_cruise=CD0_cruise,
                FuelFrac_cruise=FuelFrac_cruise,
                ThrustReduction_cruise=ThrustReduction_cruise
            )

            T_req = TW_req_SLS * W0

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

# runs cruise constraint
FuelFrac_cruise = TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction

T_cruise, W0_cruise = outer_loop_cruise_constraint(
    wing_area_grid=WingAreaGrid,
    TOGW_guess_init=TOGW_Guess,
    T_total_guess_init=TotalThrustInitialGuess,
    CD0_cruise=0.00696,   
    q_cruise=q_cruise,
    FuelFrac_cruise=FuelFrac_cruise,
    ThrustReduction_cruise=ThrustReduction,
    tol_T_rel=1e-3,
    max_iter_T=60,
    relax=0.4
)

def dash_TW_req(WS, q_dash, CD0_dash):
    return (q_dash * CD0_dash) / WS + (K * WS) / q_dash

def outer_loop_dash_constraint(
    wing_area_grid=WingAreaGrid,
    TOGW_guess_init=30000,
    T_total_guess_init=20000,
    CD0_dash=0.00696,       
    q_dash=1000.0,
    tol_T_rel=1e-3,
    max_iter_T=60,
    relax=0.4
):

    T_converged = []
    W0_converged = []

    for S_wing in wing_area_grid:
        T_total = T_total_guess_init

        for iter_outer in range(max_iter_T):
            T0_per_engine = T_total / NumberOfEngines

            W0, converged_inner, it_inner, _ = Weight_Inner_Loop(
                TOGW_Guess       = TOGW_guess_init,
                WingArea         = S_wing,
                HorizTailArea    = HorizTailArea,
                VertTailArea     = VertTailArea,
                WetAreaFuse      = WetAreaFuse,
                NumberOfEngines  = NumberOfEngines,
                WeightCrew       = WeightCrew,
                WeightPayload    = WeightPayload,
                T_0              = T0_per_engine
            )

            if not converged_inner:
                print(f"[inner loop did not converge at S={S_wing}]")
                break

            WS = W0 / S_wing
            TW_req = dash_TW_req(WS, q_dash, CD0_dash)
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

# runs dash constraint
T_dash, W0_dash = outer_loop_dash_constraint(
    wing_area_grid=WingAreaGrid,
    TOGW_guess_init=TOGW_Guess,
    T_total_guess_init=TotalThrustInitialGuess,
    CD0_dash=CD0_dash,
    q_dash=q_dash,
    tol_T_rel=1e-3,
    max_iter_T=60,
    relax=0.4
)
    
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

# Climb constraint
plt.plot(WingAreaGrid, T_climb,
         linestyle='-', linewidth=2.0, color='C0',
         label='Climb (ROC + lapse adjustment)')

# Landing constraint
plt.plot(S_converged_landing, T_grid,
         linestyle='--', linewidth=2.0, color='C3',
         label=f'Landing (max speed = {maxLandSpeed:.1f} ft/s)')

# Launch / Takeoff constraint
plt.plot(S_min_launch, T_for_plot,
         linestyle=':', linewidth=2.0, color='C4',
         label='Launch / Takeoff (min S for WS limit)')

# Sustained Turn (8 deg/s)
plt.plot(WingAreaGridManeuver, T_maneuver,
         linestyle='-.', linewidth=2.0, color='C1',
         label='Sustained Turn (8 deg/s min @ 20k ft)')

# Sustained Turn (10 deg/s)
plt.plot(WingAreaGridManeuver, T_maneuver_10,
         linestyle=':', linewidth=2.0, color='C2',
         label='Sustained Turn (10 deg/s desired @ 20k ft)')

# Cruise constraint
plt.plot(WingAreaGrid, T_cruise,
         linestyle='-', linewidth=2.0, color='C5',
         label='Cruise (30k ft)')

# Dash constraint
plt.plot(WingAreaGrid, T_dash,
         linestyle='-', linewidth=2.0, color='C6',
         label='Dash (Mach 2 @ 30k ft)')

# Example aircraft (keep stars only)
plt.scatter(S_F22, T_F22, color='k', marker='*', s=180,
            label=f'F-22 Raptor (S={S_F22:.0f} ft², T={T_F22:.0f} lbf)')
plt.scatter(S_F35, T_F35, color='m', marker='*', s=180,
            label=f'F-35 Lightning II (S={S_F35:.1f} ft², T={T_F35:.0f} lbf)')
plt.scatter(S_DRM, T_DRM, color='r', marker='*', s=180,
            label=f'Rafale M (S={S_DRM:.1f} ft², T={T_DRM:.0f} lbf)')
plt.scatter(S_F18, T_F18, color='b', marker='*', s=180,
            label=f'F-18 Super Hornet (S={S_F18:.0f} ft², T={T_F18:.0f} lbf)')

# Our chosen design point 
plt.scatter(S_Design, T_Design, color='g', marker='*', s=180,
            label=f'Our Design Point (S={S_Design:.0f} ft², T={T_Design:.0f} lbf)')

# Formatting
plt.xlabel("Wing Area S  (ft²)", fontsize=13)
plt.ylabel("Total Thrust Required  (lbf)", fontsize=13)
plt.title("Aircraft Sizing Constraint Diagram", fontsize=15, fontweight='bold')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=9, loc='upper left', framealpha=0.95, edgecolor='gray',
           ncol=1, columnspacing=1.0, labelspacing=0.5)
plt.tight_layout()
plt.show()