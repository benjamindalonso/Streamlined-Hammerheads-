import numpy as np
import matplotlib.pyplot as plt
import math

S_wing = 400
S_ht = 40
S_vt = 30
S_wet_fuselage = 500
num_engines = 1
AR = 2.5

Wcrew = 400
Wpayload = 6802

L_D_max = 14
R = 2000          # nmi
c = 0.8          # 1/hr
V = 550           # nmi/hr
E = 0.333         # hr loiter

Clmax = 1.5 # Maximum coefficient of lift (this will occur right at stall - max angle of attack)
ROC = 200 # Rate of climb ft/min
V_horizontal = 500 * 1.68781  # Climb airspeed knots to ft/s
V_horizontal_min = V_horizontal * 60  # Climb airspeed in ft/min
Climb_Gradient = ROC / V_horizontal_min # Climb gradient
Climb_Cd0 = 0.01696 # Climb drag coefficient (I just used the openVSP value again)
e = 0.8 # Oswald efficiency factor (typical value for fighters)
Ks = 1.8 # Stall speed factor 
K = 1/(math.pi*e*AR) # Induced drag factor
rhoTropicalDay = 0.00219 # Air density at sea level on a tropical day in slugs/ft^3 (for launch constraint)
maxLandSpeed = 202.6 # Max landing speed in feet per second 



def Climb_Constraint(Ks, K, Climb_Cd0, Clmax, Climb_Gradient):
    Climb_Intial = ((Ks**2*Climb_Cd0)/(Clmax))+((K*((Clmax)/Ks**2))+((Climb_Gradient)))
    Climb_Constraint = Climb_Intial*((1/.8)*(1/.99)) # Adjusting for fuel fractions and thrust reduction
    return Climb_Constraint

# Landing Constraint finds wing loading for landing weights lbf,speeds ft/s, density slug/ft^3
def Landing_Constraint(GTOW,landWeight,maxLandSpeed,CLmaxLand,density):
    S_land = 2*landWeight/(density*((maxLandSpeed/1.15)**2)*CLmaxLand)
    landWingLoading = GTOW/S_land
    return landWingLoading

def calc_e_w(T0):
    Wdry=.521*T0**.9
    woil=.082*T0**.65
    wrev=.034*T0
    wcontrol=.26*T0**.5
    wstart=9.33*(Wdry/1000)**1.078
    weng=Wdry+woil+wrev+wcontrol+wstart
    return weng

def calc_empty_w(S_wing, S_ht, S_vt, S_wet_fuselage, TOGW, T0, num_engines):
    Wwing=S_wing*19
    wht=S_ht*4
    wvt=S_vt*5.3
    wfus=S_wet_fuselage*4.8
    Wlg=.033*TOGW
    Engw=calc_e_w(T0)
    Wengines=Engw*num_engines*1.3
    Wall=.17*TOGW
    Wempt=Wwing+wht+wvt+wfus+Wlg+Wengines+Wall
    return Wempt

def calc_wf(L_D_max, R, E, c, V):
    Ld=.94*L_D_max
    cruise = np.exp((-R * c) / (V * Ld))
    loiter = np.exp((-E * c) / (Ld))
    warmup   = 0.99
    taxi = .99
    takeoff  = 0.99
    climb    = 0.96
    combat   = 0.94
    descent  = 0.99
    landing  = 0.995
    goaround = 0.99
    FinalFuel= cruise*loiter*warmup*taxi*takeoff*climb*combat*descent*landing*goaround
    TotalFuel=(1-FinalFuel)*1.06
    return TotalFuel

def innerloopweight(TOGW_guess, S_wing, S_ht, S_vt, S_wet_fuselage, num_engines, Wcrew, Wpayload, T0, err=1e-6, max_iter=200):
    W0_history = []
    delta = np.inf
    it = 0

    while delta > err and it < max_iter:
        TotalFuel = calc_wf(L_D_max, R, E, c, V)
        W_fuel = TotalFuel*TOGW_guess

        Wempt=calc_empty_w(S_wing, S_ht, S_vt, S_wet_fuselage, TOGW_guess, T0, num_engines)

        W0_new = Wempt+Wcrew+Wpayload+W_fuel
        W0_history.append(W0_new)

        delta=abs(W0_new-TOGW_guess)/max(abs(W0_new), 1e-9)

        TOGW_guess = W0_new
        it += 1
        
    converged = (delta <= err)
    return TOGW_guess, converged, it, np.array(W0_history)
    
W0, conv, it, hist = innerloopweight(
    TOGW_guess=20000,
    S_wing=S_wing,
    S_ht=S_ht,
    S_vt=S_vt,
    S_wet_fuselage=S_wet_fuselage,
    num_engines=num_engines,
    Wcrew=Wcrew,
    Wpayload=Wpayload,
    T0=5000
)

plt.figure()
plt.plot(hist, marker='o')
plt.xlabel("Iteration")
plt.ylabel("Takeoff Gross Weight W0")
plt.title("TOGW Convergence")
plt.grid(True)
plt.show()

def outer_loop_thrust_for_climb_constraint(S_grid, TOGW_guss_init, T_total_guess_init, num_engines, S_ht, S_vt, S_wet_fuselage, Wcrew, Wpayload, coef_1_climb_constraint, coef_2_climb_constraint, tol_T_rel=1e-3, max_iter_T = 100, relax=1):
    T_total_converged = []
    W0_converged = []
    iter_counts = []
    T_total_history_allS = []

    for S_wing in S_grid:

        T_total = T_total_guess_init
        T_hist = []

        for k in range(max_iter_T):
            T0 = T_total / num_engines
            W_0, wconv, it_w, W0_hist = innerloopweight(TOGW_guss_init, S_wing, S_ht, S_vt, S_wet_fuselage, num_engines, Wcrew, Wpayload, T0)

            WS = W_0 / S_wing
            TW_req = coef_1_climb_constraint / WS + coef_2_climb_constraint*WS
            T_req = TW_req * W_0
            T_hist.append(T_total)

            if abs(T_req - T_total) / max(abs(T_total), 1e-9) < tol_T_rel:
                T_total = T_req
                break

            T_total = (1-relax)*T_total +relax*T_req

        T_total_converged.append(T_total)
        W0_converged.append(W_0)
        iter_counts.append(k+1)
        T_total_history_allS.append(T_hist)

    return (np.array(T_total_converged), np.array(W0_converged), np.array(iter_counts), T_total_history_allS, W_0, wconv, it_w, W0_hist)


V_fts = V * 6076.12 / 3600.0     
rho_cruise = 0.0009              
q_cruise = 0.5 * rho_cruise * V_fts**2
CD0_cruise = 0.00696             

def cruise_constraint_TW_req(WS, q_cruise, CD0_cruise, AR, e):
    k = 1.0 / (math.pi * e * AR)
    return (q_cruise * CD0_cruise) / WS + (k * WS) / q_cruise


def outer_loop_thrust_for_cruise_constraint(
    S_grid, TOGW_guss_init, T_total_guess_init,
    num_engines, S_ht, S_vt, S_wet_fuselage,
    Wcrew, Wpayload,
    q_cruise, CD0_cruise, AR, e,
    tol_T_rel=1e-3, max_iter_T=100, relax=1.0
):
    T_total_converged, W0_converged, iter_counts, T_total_history_allS = [], [], [], []

    for S_wing in S_grid:
        T_total = T_total_guess_init
        T_hist = []

        for k_it in range(max_iter_T):
            T0 = T_total / num_engines

            W_0, wconv, it_w, W0_hist = innerloopweight(
                TOGW_guss_init, S_wing, S_ht, S_vt, S_wet_fuselage,
                num_engines, Wcrew, Wpayload, T0
            )

            WS = W_0 / S_wing
            TW_req = cruise_constraint_TW_req(WS, q_cruise, CD0_cruise, AR, e)
            T_req = TW_req * W_0

            T_hist.append(T_total)

            if abs(T_req - T_total) / max(abs(T_total), 1e-9) < tol_T_rel:
                T_total = T_req
                break

            T_total = (1 - relax) * T_total + relax * T_req

        T_total_converged.append(T_total)
        W0_converged.append(W_0)
        iter_counts.append(k_it + 1)
        T_total_history_allS.append(T_hist)

    return (np.array(T_total_converged),
            np.array(W0_converged),
            np.array(iter_counts),
            T_total_history_allS)


def landing_constraint_wing_area(T_grid,TOGW_guss_init,maxLandSpeed,CLmaxLand,density):
    tolerance = 10**(-6)
    S_converged = []
    for T in T_grid:
        S_wing = 1
        delta = 2*tolerance
        while delta > tolerance:
             W_0, wconv, it_w, W0_hist = innerloopweight(TOGW_guss_init, S_wing, S_ht, S_vt, S_wet_fuselage, num_engines, Wcrew, Wpayload, T)
             fuelWeightLand = (W_0 - (calc_empty_w(S_wing,S_ht,S_vt,S_wet_fuselage,W_0,T,1) + Wcrew + Wpayload))*0.25
             landWeight =  fuelWeightLand + calc_empty_w(S_wing,S_ht,S_vt,S_wet_fuselage,W_0,T,1) + Wcrew + Wpayload
             W_Sreq = Landing_Constraint(W_0,landWeight,maxLandSpeed,CLmaxLand,density)
             Snew = W_0/W_Sreq
             delta = abs(Snew - S_wing)
             S_wing = Snew
        S_converged.append(Snew)
    return np.array(S_converged)   

S_grid = np.linspace(300, 600, 7)
T_grid = np.linspace(0,70000,10)

tconv, W0conv, iters, T_hist_allS, *_ = outer_loop_thrust_for_climb_constraint(
    S_grid=S_grid,
    TOGW_guss_init=40000,
    T_total_guess_init=20000,
    num_engines=2,
    S_ht=S_ht,
    S_vt=S_vt,
    S_wet_fuselage=S_wet_fuselage,
    Wcrew=Wcrew,
    Wpayload=Wpayload,
    coef_1_climb_constraint=0.03,
    coef_2_climb_constraint=0.0004
)

S_convereged_landing_constraint = landing_constraint_wing_area(T_grid,40000,maxLandSpeed,Clmax,rhoTropicalDay)

plt.figure()
plt.plot(S_grid, tconv, marker='o')
plt.xlabel("Wing Area S (ft²)")
plt.ylabel("Total Thrust (lb)")
plt.title("Climb Constraint: Thrust vs Wing Area")
plt.grid(True)
plt.show()

tconv_cruise, W0conv_cruise, iters_cruise, Thist_cruise = outer_loop_thrust_for_cruise_constraint(
    S_grid=S_grid,
    TOGW_guss_init=40000,
    T_total_guess_init=20000,
    num_engines=2,
    S_ht=S_ht,
    S_vt=S_vt,
    S_wet_fuselage=S_wet_fuselage,
    Wcrew=Wcrew,
    Wpayload=Wpayload,
    q_cruise=q_cruise,
    CD0_cruise=CD0_cruise,
    AR=AR,
    e=e,
    tol_T_rel=1e-3,
    max_iter_T=100,
    relax=1.0
)

plt.figure()
plt.plot(S_grid, tconv_cruise, marker='o')
plt.xlabel("Wing Area S (ft²)")
plt.ylabel("Total Thrust (lb)")
plt.title("Cruise Constraint: Thrust vs Wing Area")
plt.grid(True)
plt.show()

plt.figure()
plt.plot(S_convereged_landing_constraint,T_grid)
plt.xlabel("Wing Area S (ft²)")
plt.ylabel("Total Thrust (lb)")
plt.title("Landing Constraint: Thrust vs Wing Area")
plt.grid(True)
plt.show()
