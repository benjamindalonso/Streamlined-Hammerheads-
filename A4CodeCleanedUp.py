import numpy as np
import matplotlib.pyplot as plt
import math

# Parameters:

# Design Point
T_Design = 43000 # lbf
S_Design = 600 # ft^2
Clean_Cd0 = 0.00766

# Weight Calculations
TOGW_Guess = 45000 # lbf
wing_area = 550 # ft^2
horiz_tail_area = 77.7 # ft^2
vert_tail_area = 86.21 # ft^2
wet_area_fuse = 407 # ft^2
number_of_engines = 1
t_0 = 43000 # Initial thrust guess per engine in lbf 
MaxLD = 13 # Maximum Lift to Drag Ratio
Range = 2000 # Range in nautical miles
Endurance = 0.333 # Loiter time in hours
Ct = 0.889 # Thrust specific fuel consumption in lb/hr/lbf
CruiseSpeed = 550 # Cruise speed in nm/hr
WingAreaGrid = list(range(200, 1500, 30)) # Grid of wing areas to evaluate in ft^2
ThrustGrid = list(range(0,100000,100)) # Grid of thrust values to evaluate in lbf

# Climb Constraint 
Ks = 1.8 # Stall speed factor
e = 0.8 # Oswald efficiency factor (Typical value for fighter))
AR = 2.227 # Aspect Ratio from OpenVSP
K =  1/(math.pi*e*AR) # Induced drag factor
Climb_Cd0 = Clean_Cd0 + 0.015 # Adjust coefficient of drag for takeoff flaps
Takeoff_Clmax = 2.0 # Slide 11 Preliminary sizing lecture part 2
ROC = 200 # Rate of _ ft/min
V_horizontal = 500 * 1.68781  # Climb airspeed knots to ft/s
V_horizontal_min = V_horizontal * 60  # Climb airspeed in ft/min
Climb_Gradient = ROC / V_horizontal_min # Climb gradient
rhoTropicalDay = 0.00219 # Air density at sea level on a tropical day in slugs/ft^3 (for landing constraint)

# Landing Constraint
maxLandSpeed = 202.6 # Max landing speed in feet per second 
Landing_Clmax = 2.3

# Launch Constraint
Vend = 135 # Catapult end speed in knots with a 67,000 GTOW and a 210 CSV setting on the catapult 
Vwod = 0 # Wind speed over the deck in knots (Assumed 0 for worst case scenario)
Vthrust = 10 # Velocity added by engine thrust during catapult launch (Assumed to be 10 knots per Raymer page 136)
Vendfps = 1.6878 * Vend # This converts knots to feet per second
Vwodfps = 1.6878 * Vwod # This converts knots to feet per second
Vthrustfps = 1.6878 * Vthrust # This converts knots to feet per second
Takeoff_Clmax = 2.0 # Clmax at takeoff per slide 11 of preliminary sizing part 2

# Stall Constraint
Vstall = 120 # Target stall speed in knots
Vstallfps = 1.6878 * Vstall # Stall speed in feet per second
Stall_Clmax = 2.6 # Maximum cl based on slide 11 of pleminiary sizing part 2

# Maneuver Constraint
TurnRate = 0.139626 # in Rad/s based on rfp minimum of 8 deg/s
TurnRate_10deg = 0.174533 # in Rad/s based on rfp preference of 10 deg/s
g_20k = 32.113 # ft/s^2 (gravitational acceleration)
Vturn = 500 # Velocity during the maneuver in knots
Vturnfps = 1.6878 * Vturn # Velocity during the maneuver in feet per second
rho_20k = 0.001267 # slug/ft^3 at 20k feet per rfp maneuver altitude
rho_30k = 0.000891 # slug/ft^3 at 30k feet 
MidMissionFuelFraction =  0.83 # Fuel fraction halfway through cruise portion of mission (This is based on the rfp requirements for maneuvering being done at "mid mission weight")
TakeoffFuelFraction = 0.99 # Fuel fraction from assignment2code
ClimbFuelFraction = 0.96 # Fuel fraction from assignment2code
ThrustReduction =  0.8 # Thrust reduction factor at cruise (due to altitude and speed)

# Cruise Constraint
Cruisefps = CruiseSpeed * 1.6878 # Cruise speed in feet per second
q_cruise = 0.5 * rho_30k * (Cruisefps**2) 
CruiseFuelFraction = TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction
CruiseThrustReduction = 0.8 # Thrust reduction factor at cruise (due to altitude and speed)

# Dash Constraint
DashSpeed = 1150 # Dash speed in knots (Mach 2 at 30k ft)
DashSpeedfps = DashSpeed * 1.6878 # Dash speed in feet per second
q_dash = 0.5 * rho_30k * (DashSpeedfps**2)



# Calculates weight of the crew
NumberOfPilots = 1 # Number of pilots
WeightOfPilot = 250 # Weight of pilot in lbf
WeightCrew = NumberOfPilots * WeightOfPilot

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
def Empty_Weight_Calculation(WingArea, HorizTailArea, VertTailArea, WetAreaFuse, TOGW_Guess, T_0 , NumberOfEngines):
    WingWeight = WingArea * 9
    HorizTailWeight = HorizTailArea * 4
    VertTailWeight = VertTailArea * 5.3
    FuseWeight = WetAreaFuse * 4.8
    LandingGearWeight = 0.045 * TOGW_Guess
    EWeight = Engine_Weight_Calculation(T_0)
    EngineWeight = EWeight * NumberOfEngines * 1.3
    AllElse = 0.17 * TOGW_Guess
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

# Weight Inner Loop Calculation Definition 
def Weight_Inner_Loop(TOGW_Guess, WingArea, HorizTailArea, VertTailArea, WetAreaFuse, NumberOfEngines, WeightCrew, WeightPayload, T_0, err=1e-6, max_iter=1000):
    W0_history = []
    delta = np.inf
    it = 0
    while delta > err and it < max_iter:
        # 1) Fuel Fraction
        Wf_W0 = Weight_Fraction_Calculation(MaxLD, Range, Endurance, Ct, CruiseSpeed)
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

W0, conv, it, hist = Weight_Inner_Loop(TOGW_Guess, wing_area,horiz_tail_area, vert_tail_area,wet_area_fuse, number_of_engines, WeightCrew, WeightPayload, t_0)


# Plots weight convergence
print('Converged TOGW =', W0, 'lbf')
plt.figure()
plt.plot(hist, marker='o')
plt.xlabel("Iteration")
plt.ylabel("Takeoff Gross Weight W0 (lbf)")
plt.title("TOGW Convergence")
plt.grid(True)
plt.show()

#=======================================================================================

# Climb Loop
def outer_loop_climb_constraint(WingAreaGrid, TOGW_Guess, t_0, tol_T_rel=1e-3, max_iter_T=60, relax=0.4):
   
    T_Converged_Climb = [] # Allocates space for converged thrust values for climb constraint
    TOGW_Converged_Climb = [] # Allocates space for converged TOGW values for climb constraint
   
   # Calculates the required thrust-to-weight ratio for climb constraint (this is a constant)
    climb_initial = (Ks**2 * Climb_Cd0 / Takeoff_Clmax) + (K * (Takeoff_Clmax / Ks**2)) + Climb_Gradient
    TW_req = climb_initial * (1 / 0.8) * (1 / 0.99)
    
    # Loop
    for S_wing in WingAreaGrid:

        T_total = t_0 * number_of_engines # Current total thrust for iteration 
        
        for iter_outer in range(max_iter_T):
            T0_per_engine = T_total / number_of_engines
            
            W0, converged_inner, it_inner,_ = Weight_Inner_Loop(
                TOGW_Guess,
                S_wing,
                horiz_tail_area,
                vert_tail_area,
                wet_area_fuse,
                number_of_engines,
                WeightCrew,
                WeightPayload,
                T0_per_engine
            )
            
            # Check if this iteration has enough thrust to satisfy the constraint
            WS = W0 / S_wing
            T_req = TW_req * W0
            rel_error = abs(T_req - T_total) / max(abs(T_total), 1e-6)
            
            if rel_error < tol_T_rel:
                T_Converged_Climb.append(T_req)
                TOGW_Converged_Climb.append(W0)
                break
            
            T_total = (1 - relax) * T_total + relax * T_req
        
        else:
            T_Converged_Climb.append(T_req)
            TOGW_Converged_Climb.append(W0)
    
    return np.array(T_Converged_Climb), np.array(TOGW_Converged_Climb)

# Runs climb constraint
T_Converged_Climb, TOGW_Converged_Climb = outer_loop_climb_constraint(WingAreaGrid, TOGW_Guess, t_0, tol_T_rel=1e-3, max_iter_T=60, relax=0.4)


# Landing Loop
def outer_loop_landing_constraint(ThrustGrid, TOGW_Guess, maxLandSpeed, Landing_Clmax, RhoTropicalDay):
    
    tolerance = 10**(-6)
    maxLandSpeed = maxLandSpeed + 25.3 # Adds 15 kts wind over deck in ft/s
    S_Converged_Landing = []
    
    for T in ThrustGrid:
        S_wing = 1 # Starting point for wing area convergence in ft^2
        delta = 2*tolerance # Initilize delta to be larger than tolerance 
        while delta > tolerance:
             # Calculate new TOGW based on current wing area 
             W0, wconv, it_w, W0_hist = Weight_Inner_Loop(TOGW_Guess, S_wing, horiz_tail_area, vert_tail_area, wet_area_fuse, number_of_engines, WeightCrew, WeightPayload, T)
             
             # Calculate "no fuel" Weight based on current wing area and new TOGW
             EmptyWeight = Empty_Weight_Calculation(S_wing, horiz_tail_area, vert_tail_area, wet_area_fuse, W0, T, number_of_engines) + WeightCrew + WeightPayload
             
             # Calculate weight of fuel remaining on landing (25% of maximum fuel weight)
             fuelWeightLand = (W0 - EmptyWeight)*0.25

             # Calculates landing weight based on empty weight and fuel remaining on landing 
             landWeight =  fuelWeightLand + EmptyWeight

             # Calculate required wing area based on landing weight and landing speed constraint
             S_req = 2*landWeight/(RhoTropicalDay*((maxLandSpeed/1.15)**2)*Landing_Clmax)
             delta = abs(S_req - S_wing)
             S_wing = S_req
        S_Converged_Landing.append(S_wing)
    return np.array(S_Converged_Landing)  

# Runs landing constraint
S_Converged_Landing = outer_loop_landing_constraint(ThrustGrid, TOGW_Guess, maxLandSpeed, Landing_Clmax, rhoTropicalDay)

# Launch/stall loop
# Calculate the two maximum wing loadings (Launch is from Raymer page 136)
WS_max_launch = 0.5 * rhoTropicalDay * ((Vendfps + Vwodfps + Vthrustfps)**2) * Takeoff_Clmax / 1.21
WS_max_stall  = 0.5 * rhoTropicalDay * (Vstallfps ** 2) * Stall_Clmax

S_min_launch = []
T_for_plot = []
S_min_stall = []
T_for_plotstall = []

for T_total in ThrustGrid:

    T0 = T_total / number_of_engines
    
    # Get converged TOGW for each thrust level
    W0, converged, _, _ = Weight_Inner_Loop(
        TOGW_Guess,
        wing_area,         
        horiz_tail_area,
        vert_tail_area,
        wet_area_fuse,
        number_of_engines,
        WeightCrew,
        WeightPayload,
        T0,
        err=1e-5,
        max_iter=1000
    )
    
    # Minimum wing area required to meet each constraint
    S_req_launch = W0 / WS_max_launch 
    S_req_stall  = W0 / WS_max_stall  
    
   # Build arrays for plotting
    S_min_launch.append(S_req_launch)
    S_min_stall.append(S_req_stall)
    T_for_plot.append(T_total)
    T_for_plotstall.append(T_total)

    
 

# Maneuver Constraint Equation Definition
def Sustained_Turn_Constraint(TurnRate, g, Vturnfps, Clean_Cd0, K, Wing_Loading, rho_30k, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction):
    n = math.sqrt(((TurnRate * Vturnfps) / g)**2 + 1)
    q = 0.5 * rho_30k * (Vturnfps**2)
    Wing_LoadingTurn = Wing_Loading * TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction 
    
    # Required T/W at maneuver alt
    TW_turn = (q * Clean_Cd0 / Wing_LoadingTurn) + (K * (n**2) * Wing_LoadingTurn / q)
    
    # Correct back to Sea Level Static T/W (T0/W0)
    TW_SLS = TW_turn * (TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction / ThrustReduction)
    return TW_SLS

# Maneuver Loop
def outer_loop_maneuver_constraint(
        WingAreaGrid,
        TOGW_Guess,
        ThrustPerEngine,      
        TurnRate_input,
        tol_T_rel=1e-3,       
        max_iter_T=60,
        relax=0.4):           

    T_Converged_Maneuver = []

    for WingArea in WingAreaGrid:

        T_total = ThrustPerEngine * number_of_engines   # starting guess

        for k in range(max_iter_T):

            T0_per_engine = T_total / number_of_engines

            # Inner weight loop
            Final_TOGW, Converged, Iterations, _ = Weight_Inner_Loop(
                TOGW_Guess,
                WingArea,
                horiz_tail_area,      
                vert_tail_area,
                wet_area_fuse,
                number_of_engines,
                WeightCrew,
                WeightPayload,
                T0_per_engine         
            )

            WingLoading = Final_TOGW / WingArea

            # Get required T/W from turn function
            TW_maneuver = Sustained_Turn_Constraint(
                TurnRate_input, g_20k, Vturnfps, Clean_Cd0, K,
                WingLoading, rho_20k,
                MidMissionFuelFraction, TakeoffFuelFraction,
                ClimbFuelFraction, ThrustReduction
            )
            T_required = TW_maneuver * Final_TOGW
            if T_required > 100000:
                T_required = 100000 # Cap required thrust to avoid numerical issues in convergence

            rel_error = abs(T_required - T_total) / max(abs(T_total), 1e-9)

            if rel_error < tol_T_rel:
                break

            # Relaxation (same as your climb constraint)
            T_total = (1 - relax) * T_total + relax * T_required

        # After inner iterations, store the result
        T_Converged_Maneuver.append(T_total)

    return np.array(T_Converged_Maneuver)

# Run Maneuver Constraint Loop
T_Converged_Maneuver = outer_loop_maneuver_constraint(WingAreaGrid, TOGW_Guess, t_0, TurnRate, tol_T_rel=1e-3, max_iter_T=60, relax=0.4)
    
# Run Maneuver Constraint Loop For 10 deg/s
T_Converged_Maneuver10 = outer_loop_maneuver_constraint(WingAreaGrid, TOGW_Guess, t_0, TurnRate_10deg, tol_T_rel=1e-3, max_iter_T=60, relax=0.4)
    
def cruise_TW_req_SLS(WingLoading_W0S, q_cruise, CD0_cruise,
                      FuelFrac_cruise, ThrustReduction_cruise):
    WS_cruise = WingLoading_W0S * FuelFrac_cruise
    TW_alt = (q_cruise * CD0_cruise) / WS_cruise + (K * WS_cruise) / q_cruise
    TW_SLS = TW_alt * (FuelFrac_cruise / ThrustReduction_cruise)
    return TW_SLS

# Cruise Constraint Loop
def outer_loop_cruise_constraint(
   WingAreaGrid,
   TOGW_Guess,
   t_0,
   Clean_Cd0,    
    q_cruise,
    FuelFrac_cruise, 
    ThrustReduction_cruise,
    tol_T_rel=1e-3,
    max_iter_T=60,
    relax=0.4
):
    T_Converged_Cruise = []

    for S_wing in WingAreaGrid:

        for iter_outer in range(max_iter_T):

            T_total = t_0 * number_of_engines  # Reset thrust for this wing area

            W0, converged_inner, it_inner, _ = Weight_Inner_Loop(
                TOGW_Guess,
                S_wing,
                horiz_tail_area,
                vert_tail_area,
                wet_area_fuse,
                number_of_engines,
                WeightCrew,
                WeightPayload,
                t_0
            )

            WS_W0 = W0 / S_wing  # W0/S

            TW_req_SLS = cruise_TW_req_SLS(
                WS_W0,
                q_cruise,
                Clean_Cd0,
                CruiseFuelFraction,
                CruiseThrustReduction
            )

            T_req = TW_req_SLS * W0

            rel_error = abs(T_req - T_total) / max(abs(T_total), 1e-6)
            if rel_error < tol_T_rel:
                T_Converged_Cruise.append(T_req)
                break

            t_0 = (1 - relax) * t_0 + relax * T_req

        else:
            T_Converged_Cruise.append(T_req)

    return np.array(T_Converged_Cruise)

# runs cruise constraint

T_Converged_Cruise = outer_loop_cruise_constraint(
    WingAreaGrid,
    TOGW_Guess,
    t_0,
    Clean_Cd0,   
    q_cruise,
    CruiseFuelFraction,
    CruiseThrustReduction,
    tol_T_rel=1e-3,
    max_iter_T=60,
    relax=0.4
)

# Runs Dash Constraint Using Cruise Loop With Different Input Data

T_Converged_Dash = outer_loop_cruise_constraint(
    WingAreaGrid,
    TOGW_Guess,
    t_0,
    Clean_Cd0,   
    q_dash,
    CruiseFuelFraction,
    CruiseThrustReduction,
    tol_T_rel=1e-3,
    max_iter_T=60,
    relax=0.4
)




# Fighter T vs S comparison data

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



# =============================================================================
#   FINAL CONSTRAINT DIAGRAM + SUMMARY OUTPUT
# ============================================================================


# ─── Combined Constraint Diagram ────────────────────────────────────────────
plt.figure(figsize=(11, 7))

# Climb constraint
plt.plot(WingAreaGrid, T_Converged_Climb,
         linestyle='-', linewidth=2.0, color='C0',
         label='Climb (ROC + lapse adjustment)')

# Landing constraint
plt.plot(S_Converged_Landing, ThrustGrid,
         linestyle='--', linewidth=2.0, color='C3',
         label=f'Landing (max speed = {maxLandSpeed:.1f} ft/s)')

# Launch / Takeoff constraint
plt.plot(S_min_launch, T_for_plot,
         linestyle=':', linewidth=2.0, color='C4',
         label='Launch / Takeoff (min S for WS limit)')

# Stall Constraint
plt.plot(S_min_stall, T_for_plotstall,
         linestyle=':', linewidth=2.0, color='C8',
         label='Stall (min S for WS limit)')

# Sustained Turn (8 deg/s)
plt.plot(WingAreaGrid, T_Converged_Maneuver,
         linestyle='-.', linewidth=2.0, color='C1',
         label='Sustained Turn (8 deg/s min @ 30k ft)')

# Sustained Turn (10 deg/s)
plt.plot(WingAreaGrid, T_Converged_Maneuver10,
         linestyle=':', linewidth=2.0, color='C2',
         label='Sustained Turn (10 deg/s desired @ 20k ft)')

# Cruise constraint
plt.plot(WingAreaGrid, T_Converged_Cruise,
         linestyle='-', linewidth=2.0, color='C5',
         label='Cruise (30k ft)')

# Dash constraint
plt.plot(WingAreaGrid, T_Converged_Dash,
         linestyle='-', linewidth=2.0, color='C6',
         label='Dash (Mach 2 @ 30k ft)')




# Example aircraft (keep stars only)
plt.scatter(S_F22, T_F22, color='k', marker='o', s=50,
            label=f'F-22 Raptor (S={S_F22:.0f} ft², T={T_F22:.0f} lbf)')
plt.scatter(S_F35, T_F35, color='m', marker='o', s=50,
            label=f'F-35 Lightning II (S={S_F35:.1f} ft², T={T_F35:.0f} lbf)')
plt.scatter(S_DRM, T_DRM, color='g', marker='o', s=50,
            label=f'Rafale M (S={S_DRM:.1f} ft², T={T_DRM:.0f} lbf)')
plt.scatter(S_F18, T_F18, color='b', marker='o', s=50,
            label=f'F-18 Super Hornet (S={S_F18:.0f} ft², T={T_F18:.0f} lbf)')

# Our chosen design point 
plt.scatter(S_Design, T_Design, color='r', marker='*', s=300,
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