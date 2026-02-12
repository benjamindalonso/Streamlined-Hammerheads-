# This script calculates the launch constraint on a T vs S plot

import numpy as np
import matplotlib.pyplot as plt
import math

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

# Maneuvering Constraint (sustained and instant)
def Sustained_Turn_Constraint(TurnRate, g, Vturn, Cd0, K, Wing_Loading, rho, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction):
    n = math.sqrt(((TurnRate * Vturn) / g)**2 + 1)
    q = 0.5 * rho * (Vturn**2)
    Wing_LoadingTurn = Wing_Loading * TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction 
    
    # Required T/W at maneuver alt
    TW_turn = (q * Cd0 / Wing_LoadingTurn) + (K * (n**2) * Wing_LoadingTurn / q)
    
    # Correct back to Sea Level Static T/W (T0/W0)
    TW_SLS = TW_turn * (TakeoffFuelFraction * ClimbFuelFraction * MidMissionFuelFraction / ThrustReduction)
    return TW_SLS

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

# Fixed parameters
MaxLD = 12
Range = 2000                  # nmi
Endurance = 0.333            # hr
Ct = 0.889                   # lb/(lbf hr)
V = 550                      # knots
HorizTailArea = 100
VertTailArea = 81
WetAreaFuse = 180
NumberOfEngines = 1  # Example number of engines
TurnRate = 0.1745 # in Rad/s based on rfp preference of 10 deg/s
g = 32.174 # ft/s^2 (gravitational acceleration)
Vturn = 750 # Velocity during the maneuver in feet per second
Cd0Turn = 0.00696 # Zero lift drag coefficient during the turn (I just used the openVSP value again)
e = 0.8
AR = 2.5
K =  1/(math.pi*e*AR) # Induced drag factor
rhoTurn = 0.001267 # Slugs per ft^3 (Density at 20,000 ft - maneuvering altitude per rfp)
MidMissionFuelFraction =  0.906 # Fuel fraction half way through cruise portion of mission (This is based on the rfp requirements for maneuvering being done at "mid mission weight")
TakeoffFuelFraction = 0.99 # Fuel fraction from assignment2code
ClimbFuelFraction = 0.96 # Fuel fraction from assignment2code
ThrustReduction =  0.8 # Thrust reduction factor at cruise (due to altitude and speed)
WingArea = 610


T_0 = 43000  # Example value for thrust per engine
TotalThrustInitialGuess = T_0 * NumberOfEngines
TOGW_Guess = 30000  # Initial guess for Takeoff Gross Weight in pounds
WingAreaGrid = list(range(200, 1000, 5))  # range of wing areas to analyze


# Outer Loop
def Thrust_Outer_Loop(WingAreaGrid, TOGW_Guess, TotalThrustInitialGuess, NumberOfEngines, HorizTailArea, VertTailArea, WetAreaFuse, WeightCrew, WeightPayload, tol_T_rel=1e-3, max_iter_T=100, relax=0.3):
    T_total_converged = []
    W0_converged = []
    iter_counts = []
    T_total_history_allS = []  

    for WingArea in WingAreaGrid:
     T_total = TotalThrustInitialGuess
     T_hist = []
     for k in range(max_iter_T):
         T_0 = T_total/NumberOfEngines
         Final_TOGW, Converged, Iterations, W0_history = Weight_Inner_Loop( TOGW_Guess, WingArea, HorizTailArea, VertTailArea, WetAreaFuse, NumberOfEngines, WeightCrew, WeightPayload, T_0)
         WingLoading = Final_TOGW / WingArea
         # Maneuver Constraint
         ThrustToWeightManeuver = Sustained_Turn_Constraint(TurnRate, g, Vturn, Cd0Turn, K, WingLoading, rhoTurn, MidMissionFuelFraction, TakeoffFuelFraction, ClimbFuelFraction, ThrustReduction)
         ThrustRequiredManeuver = ThrustToWeightManeuver * Final_TOGW
         T_hist.append(T_total)
         if abs(ThrustRequiredManeuver - T_total) / max(abs(T_total), 1e-9) < tol_T_rel:
                T_total = ThrustRequiredManeuver
                break
         T_total = (1 - relax) * T_total + relax * ThrustRequiredManeuver
     T_total_converged.append(T_total)
     W0_converged.append(Final_TOGW)
     iter_counts.append(k+1)
     T_total_history_allS.append(np.array(T_hist))
    return (np.array(T_total_converged), np.array(W0_converged), np.array(iter_counts), T_total_history_allS)

T_total_curve, W0_curve, n_iter_T, T_hist_allS = Thrust_Outer_Loop(WingAreaGrid, TOGW_Guess, TotalThrustInitialGuess, NumberOfEngines, HorizTailArea, VertTailArea, WetAreaFuse, WeightCrew, WeightPayload, tol_T_rel=1e-3, max_iter_T=100, relax=1.0)

plt.figure(figsize=(16,9))
plt.title('Converged T vs S for Maneuver Constraint')
plt.xlabel("Wing Area S (ft^2)")
plt.ylabel("Total Thrust T (lbf)")
plt.plot(WingAreaGrid, T_total_curve, label='Converged T for Maneuver Constraint', marker='o')
plt.legend(loc='best')
plt.grid()
plt.show()