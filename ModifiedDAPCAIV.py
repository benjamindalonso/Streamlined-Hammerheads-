# This script runs a modified DAPCA IV preliminary cost analysis for an aircraft design.

# Equations
def engineering_hours(W_Empty, V, Q):
    """Engineering hours HE (imperial/fps units)"""
    return (4.86 * (W_Empty ** 0.777) * (V ** 0.894) * (Q ** 0.163)) * HoursFudgeFactor * StealthFudgeFactor


def tooling_hours(W_Empty, V, Q):
    """Tooling hours HT"""
    return (5.99 * (W_Empty ** 0.777) * (V ** 0.696) * (Q ** 0.263)) * HoursFudgeFactor * StealthFudgeFactor


def manufacturing_hours(W_Empty, V, Q):
    """Manufacturing hours HM"""
    return (7.37 * (W_Empty ** 0.82) * (V ** 0.484) * (Q ** 0.641)) * HoursFudgeFactor * StealthFudgeFactor


def qc_hours(HM):
    """QC hours HQ non-cargo (0.133 factor)"""
    return (0.133 * HM) * HoursFudgeFactor * StealthFudgeFactor


def development_support_cost(W_Empty, V):
    """Development support cost CD (Dollars)"""
    return 91.3 * (W_Empty ** 0.630) * (V ** 1.3) * StealthFudgeFactor * InflationFudgeFactor


def flight_test_cost(W_Empty, V, FTA):
    """Flight test cost CF (Dollars)"""
    return 2498 * (W_Empty ** 0.325) * (V ** 0.822) * (FTA ** 1.21) * StealthFudgeFactor * InflationFudgeFactor


def manufacturing_materials_cost(W_Empty, V, Q):
    """Manufacturing materials cost CM (Dolalrs)"""
    return 22.1 * (W_Empty ** 0.921) * (V ** 0.621) * (Q ** 0.799) * StealthFudgeFactor * InflationFudgeFactor


def engine_production_cost(T_Max, M_Max, T_TurbInlet):
    """Engine production cost per enigne (Dollars)"""
    inner = (
        0.043 * T_Max +
        243.25 * M_Max +
        0.969 * T_TurbInlet -
        2228
    )
    return (3112 * inner) * StealthFudgeFactor * InflationFudgeFactor


def total_rdte_and_flyaway(
    HE, HT, HM, HQ,
    RE, RT, RM, RQ,
    CD, CF, CM, Ceng_per_engine, N_Eng
):
    """Total RDT&E + flyaway cost (Dollars)"""
    labor_cost = (
        HE * RE +
        HT * RT +
        HM * RM +
        HQ * RQ
    )
    
    non_labor = (
        CD +
        CF +
        CM +
        (Ceng_per_engine * N_Eng)
    )
    
    return labor_cost + non_labor



# Input Parameters
W_Empty = 67822            # Empty Weight lbs
V = 900                  # Max Velocity kts (nm/hr)
Q = 500                     # Lesser of production quantity or number to be produced in five years
FTA = 2                     # Number of flight-test aircraft (usually 2-6)
N_Eng = 500                 # Total production quantity times number of engines per aircraft 
T_Max = 43000               # Engine maximum thrust lbs
M_Max = 1.6               # engine maximum mach number
T_TurbInlet = 4059.67       # turbine inlet temperature R
HoursFudgeFactor = 1.2      # Fudge factor to account for complex materials. Aluminum = 1, Graphite-Epoxy = 1.1-1.8, Fiberglass = 1.1-1.2, Steel = 1.5-2, Titanium = 1.1-1.8
InflationFudgeFactor = 2    # Assume 1 dollar in 2012 is worth 2 dollars in 2035
StealthFudgeFactor = 1.2    # Fudge factor to account for stealth requirements in manufacturing
InflationFudgeFactor = 2    # Assume 1 dollar in 2012 is worth 2 dollars in 2035

# Hourly wrap rates (include salaries, benefits, overhead, and administrative costs)
year = 2035                    # Target year of opperation
RT   = 2.883 * year - 5666     # Tooling
RE   = 2.576 * year - 5058     # Engineering
RQC  = 2.60  * year - 5112     # QC
RM   = 2.316 * year - 4552     # Manufacturing

# Do calculations
HE  = engineering_hours(W_Empty, V, Q)
HT  = tooling_hours(W_Empty, V, Q)
HM  = manufacturing_hours(W_Empty, V, Q)
HQ  = qc_hours(HM)
CD  = development_support_cost(W_Empty, V)
CF  = flight_test_cost(W_Empty, V, FTA)
CM  = manufacturing_materials_cost(W_Empty, V, Q)
Ceng_per_engine = engine_production_cost(T_Max, M_Max, T_TurbInlet)    # cost per engine

total_cost = total_rdte_and_flyaway(HE, HT, HM, HQ, RE, RT, RM, RQC,CD, CF, CM, Ceng_per_engine, N_Eng)

# Calculate dollar cost for each labor category
cost_engineering = HE * RE
cost_tooling     = HT * RT
cost_manufacturing = HM * RM
cost_qc          = HQ * RQC

# Estimate avionics cost to be 30 percent of the total cost
C_Avionics = (3/7)*total_cost
total_cost = total_cost + C_Avionics # Update total cost to now include avionics 

# Print results
print("\n=== DAPCA-IV COST ESTIMATE (2035 dollars) ===")
print(f"Total Program Cost (RDT&E + Flyaway for {Q} aircraft): ${total_cost:,.0f}")
print(f"Flyaway Cost per Aircraft:                    ${total_cost / Q:,.0f}")
print(f"\nBreakdown:")
print(f"Engineering hours (HE)     : {HE:,.0f}")
print(f"Tooling hours (HT)         : {HT:,.0f}")
print(f"Manufacturing hours (HM)   : {HM:,.0f}")
print(f"QC hours (HQ)              : {HQ:,.0f}")
print(f"Engineering labor cost     : ${cost_engineering:,.0f}")
print(f"Tooling labor cost         : ${cost_tooling:,.0f}")
print(f"Manufacturing labor cost   : ${cost_manufacturing:,.0f}")
print(f"QC labor cost              : ${cost_qc:,.0f}")
print(f"Total labor cost           : ${cost_engineering + cost_tooling + cost_manufacturing + cost_qc:,.0f}")
print(f"Development Support (CD)   : ${CD:,.0f}")
print(f"Flight Test (CF)           : ${CF:,.0f}")
print(f"Materials (CM)             : ${CM:,.0f}")
print(f"Engines total              : ${Ceng_per_engine * N_Eng:,.0f}")
print(f"Avionics total             : ${C_Avionics:,.0f}")