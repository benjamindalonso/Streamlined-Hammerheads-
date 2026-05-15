# This script runs a modified DAPCA IV preliminary cost analysis for an aircraft design.

# Equations based on pages 732-733 of Raymers "Aircraft Design: A Conceptual Approach" 

# Calculates engineering hours
def engineering_hours(W_Empty, V, Q):
    """Engineering hours HE (imperial/knots units)"""
    return (4.86 * (W_Empty ** 0.777) * (V ** 0.894) * (Q ** 0.163)) * HoursFudgeFactor * StealthFudgeFactor

# Calculates tooling hours
def tooling_hours(W_Empty, V, Q):
    """Tooling hours HT"""
    return (5.99 * (W_Empty ** 0.777) * (V ** 0.696) * (Q ** 0.263)) * HoursFudgeFactor * StealthFudgeFactor

# Calculates manufacturing hours
def manufacturing_hours(W_Empty, V, Q):
    """Manufacturing hours HM"""
    return (7.37 * (W_Empty ** 0.82) * (V ** 0.484) * (Q ** 0.641)) * HoursFudgeFactor * StealthFudgeFactor

# Calculates quality control hours
def qc_hours(HM):
    """QC hours HQ non-cargo (0.133 factor)"""
    return (0.133 * HM) * HoursFudgeFactor * StealthFudgeFactor

# Calculates development support cost
def development_support_cost(W_Empty, V):
    """Development support cost CD (Dollars)"""
    return 91.3 * (W_Empty ** 0.630) * (V ** 1.3) * StealthFudgeFactor * InflationFudgeFactor

# Calculates flight test cost
def flight_test_cost(W_Empty, V, FTA):
    """Flight test cost CF (Dollars)"""
    return 2498 * (W_Empty ** 0.325) * (V ** 0.822) * (FTA ** 1.21) * StealthFudgeFactor * InflationFudgeFactor

# Calculates manufacturing materials cost
def manufacturing_materials_cost(W_Empty, V, Q):
    """Manufacturing materials cost CM (Dolalrs)"""
    return 22.1 * (W_Empty ** 0.921) * (V ** 0.621) * (Q ** 0.799) * StealthFudgeFactor * InflationFudgeFactor

# Calculates engine production cost per engine
def engine_production_cost(T_Max, M_Max, T_TurbInlet):
    """Engine production cost per enigne (Dollars)"""
    inner = (
        0.043 * T_Max +
        243.25 * M_Max +
        0.969 * T_TurbInlet -
        2228
    )
    return (3112 * inner) * StealthFudgeFactor * InflationFudgeFactor

# Calculates total RDT&E + flyaway cost
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
W_Empty = 31862                # Empty Weight lbs (38,199 for design B, 42970 for A))
V = 985                       # Max Velocity kts (nm/hr) (900 for design B, 1290 for A)
Q = 500                        # Lesser of production quantity or number to be produced in five years
FTA = 2                        # Number of flight-test aircraft (usually 2-6)
N_Eng = 502                  # Total production quantity times number of engines per aircraft (500 for design B)
T_Max = 43000                  # Engine maximum thrust lbs (43,000 for design B, 60,000 for A)
M_Max = 1.67                   # Engine maximum mach number
T_TurbInlet = 4059.67          # Turbine inlet temperature R
HoursFudgeFactor = 1.5         # Fudge factor to account for complex materials. Aluminum = 1, Graphite-Epoxy = 1.1-1.8, Fiberglass = 1.1-1.2, Steel = 1.5-2, Titanium = 1.1-1.8 (1.1 for design B, 1.5 for A)
InflationFudgeFactor = 1.40    # Assume 1 dollar in 2012 is worth 1.40 dollars in 2025
StealthFudgeFactor = 1.2       # Fudge factor to account for stealth requirements in manufacturing


# Hourly wrap rates (include salaries, benefits, overhead, and administrative costs)
# Based on lecture slide 20 in cost analysis presentation
year = 2025                   # Year dollars
RT   = 2.883 * year - 5666     # Tooling
RE   = 2.576 * year - 5058     # Engineering
RQC  = 2.60  * year - 5112     # QC
RM   = 2.316 * year - 4552     # Manufacturing
print(f"Hourly wrap rates (2025 dollars): RE = ${RE:,.0f}, RT = ${RT:,.0f}, RQC = ${RQC:,.0f}, RM = ${RM:,.0f}")

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

# Calculate dollar cost for each labor category (hours * hourly wrap rate)
cost_engineering = HE * RE
cost_tooling     = HT * RT
cost_manufacturing = HM * RM
cost_qc          = HQ * RQC


# Estimate avionics cost to be 30 percent of the total cost
C_Avionics = (3/7)*total_cost
total_cost = total_cost + C_Avionics # Update total cost to now include avionics 

# Seperate RDT&E cost
RDT_E_cost = (
    HE * RE +          # Engineering labor
    HT * RT +          # Tooling labor
    CD +               # Development support
    CF                 # Flight test
    + 0.2 * C_Avionics  # Assume 20% of avionics cost is RDT&E
)

# Print results
print("\n=== DAPCA-IV COST ESTIMATE (2025 dollars) ===")
print(f"Total Program Cost (RDT&E + Flyaway for {Q} aircraft): ${total_cost:,.0f}")
print(f"Flyaway Cost per Aircraft:                    ${total_cost / Q:,.0f}")
print(f"RDT&E Cost:                                   ${RDT_E_cost:,.0f}")
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


# Direct Operating Cost Estimation
FlightHoursPerYearPerAircraft = 400 # Range between 300 and 500 per Raymer
NumCrew = 1.1 * Q # Assume 1.1 crew members per aircraft per Raymer
MaintenanceManHourPerFlightHour = 12.5 # Range between 10 and 15 per Raymer

# Calculate cost of crew on duty to fly the existing aircraft
AnnualCostCrew = NumCrew * (RE * 2080) # According to Raymer

# Calculate maintenance cost per year
Ca = (total_cost / Q) - Ceng_per_engine
Ce = Ceng_per_engine
Ne = N_Eng / Q # Number of engines per aircraft
MatCostPerFH = 3.3 * (Ca / 1000000) + 14.2 + (58 * (Ce / 1000000) - 26.1) * Ne

AnnualMatCost = MatCostPerFH * FlightHoursPerYearPerAircraft * Q 

AnnualMaintenanceCost = MaintenanceManHourPerFlightHour * FlightHoursPerYearPerAircraft * Q * RM

# Here we ignore the airframe and engine depreciation due to lack of good sources

# Calculate the fuel operating costs
Vavg = 500 # knots
SFC = 0.7 # lbm/hr/lbf
Tavg = 20000 # lbf
FuelMassPerAircraftPerYear = FlightHoursPerYearPerAircraft * SFC * Tavg
TotalFuelMassPerYear = FuelMassPerAircraftPerYear * Q
TotalFuelGallonsPerYear = TotalFuelMassPerYear / 6.7 # Convert lbm to gallons
CostPerGallon = 6 # 2025 price per gallon
TotalAnnualFuelCost = TotalFuelGallonsPerYear * CostPerGallon

DOC = AnnualCostCrew + AnnualMaintenanceCost + AnnualMatCost + TotalAnnualFuelCost
print(f"\n=== DIRECT OPERATING COST ESTIMATE (2025 dollars) ===")
print(f"Annual Cost of Crew: ${AnnualCostCrew:,.0f}")
print(f"Annual Maintenance Cost: ${AnnualMaintenanceCost:,.0f}")
print(f"Annual Material Cost: ${AnnualMatCost:,.0f}")
print(f"Annual Fuel Cost: ${TotalAnnualFuelCost:,.0f}")
print(f"Total Annual Direct Operating Cost: ${DOC:,.0f}")