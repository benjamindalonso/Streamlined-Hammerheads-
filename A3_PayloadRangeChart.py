import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Aircraft Weight Data
MTOW = 54748.1   # lb  Maximum takeoff weight
OEW  = 31862.2   # lb  Operating empty weight
# Assumed from MZFW fraction of FA-18E/F (0.7438) multiplied by MTOW (54748.1 lb)
MZFW =  40732.6  # lb  Maximum zero-fuel weight (structural limit! of the aircraft)
# Maximum payload weight
MPW = MZFW - OEW       # Maximum Payload Weight  (MZFW structural limit)
 

# Maximum fuel capacity 
fuel_capacity_ft3 = 420    # ft^3
fuel_density      = 50.39  # lb/ft^3
# Calculate maximum fuel weight
MFW = fuel_capacity_ft3 * fuel_density

# Mission/Payload 
# Calculate design payload weight
DPW = 1836   # Design Payload Weight

# Reserve fuel fraction — kept on board at landing (not burned for range)
RFF = 0.15   # 15 %

print(f'  MTOW = {MTOW:>10,} lb')
print(f'  OEW  = {OEW:>10,} lb')
print(f'  Max Zero Fuel Weight = {MZFW:>10,} lb')
print(f'  Max Fuel Weight  = {MFW:>10,.0f} lb   ({fuel_capacity_ft3:,} ft^3)')
print(f'  Max Payload Weight  = {MPW:>10,} lb   (= MZFW - OEW)')
print(f'  Design Payload Weight  = {DPW:>10,} lb')
print()
print(f'  Weight budget (MTOW - OEW) = {MTOW-OEW:,.1f} lb  <- shared by payload + fuel')


# Cruise Parameters 
rho_c  = 0.000891   # slug/ft3  air density at cruise altitude
S      = 600.0      # ft2       wing reference area
CL     = 0.64180    #           cruise lift coefficient
LD     = 14.259     #           cruise L/D
CD     = 0.04501    #           cruise drag coefficient
TSFC   = 0.889      # lb/lb/hr  thrust-specific fuel consumption
ct     = TSFC / 3600# lb/lb/s   (convert for the Breguet integral)

print(f'  Altitude : [30,000] ft  ->  rho = {rho_c:.6f} slug/ft3')
print(f'  S        : {S:,.0f} ft2')
print(f'  CL = {CL},  CD = {CD:.5f},  L/D = {LD}')
print(f'  TSFC     : {TSFC} lb/lb/hr  =  {ct:.4e} lb/lb/s')


# Breguet range function
def breguet(W0, W1):
    """
    Jet Breguet range (constant altitude, constant CL).
    R = 2*sqrt(2/rho*S) * (1/ct) * (CL^0.5/CD) * (sqrt(W0) - sqrt(W1))
    W0, W1 in lb  ->  returns range in nautical miles
    """
    # FT to NM conversion factor
    FT_PER_NM = 6076.12  # ft / nautical mile

    # Compute range in feet using the Breguet formula
    R_ft = (2 * np.sqrt(2 / (rho_c * S))
              * (1 / ct)
              * (CL**0.5 / CD)
              * (W0**0.5 - W1**0.5))
    range_nm = R_ft / FT_PER_NM # convert ft to NM

    return range_nm


# Point Calculations

# Point A, Zero Range, Maximum Payload
payload_A = MPW
range_A   = 0

print(f'\\nPoint A')
print(f'  Payload = {payload_A:,} lb  (= MPW)')
print(f'  Range   = {range_A} nm')


# Point B, Harmonic Range, Maximum Payload at MTOW
FW_B      = MTOW - OEW - MPW      # fuel that fits after max payload, at MTOW
W_res_B   = RFF * FW_B            # reserve (not burned for range)
payload_B = MPW
W0_B      = MTOW
W1_B      = OEW + MPW + W_res_B
range_B   = breguet(W0_B, W1_B)

print(f'\\nPoint B')
print(f'  Fuel loaded  = MTOW - OEW - MPW = {MTOW} - {OEW} - {MPW} = {FW_B:,.1f} lb')
print(f'  Reserve      = {RFF:.0%} x {FW_B:,.1f} = {W_res_B:,.0f} lb')
print(f'  W0 = MTOW                       = {W0_B:,} lb')
print(f'  W1 = OEW + MPW + Wres           = {OEW} + {MPW} + {W_res_B:.0f} = {W1_B:,.0f} lb')
print(f'  Payload = {payload_B:,} lb')
print(f'  Range   = {range_B:,.0f} nm')


# Point C, Maximum-Fuel Range
FW_C      = MFW
payload_C = MTOW - OEW - MFW 
W_res_C   = RFF * FW_C
W0_C      = MTOW
W1_C      = OEW + payload_C + W_res_C
range_C   = breguet(W0_C, W1_C)

print(f'\\nPoint C')
print(f'  Fuel loaded  = MFW             = {FW_C:,.0f} lb')
print(f'  Payload_C    = MTOW-OEW-MFW    = {MTOW}-{OEW}-{FW_C:.0f} = {payload_C:,.0f} lb  (reduced!)')
print(f'  Reserve      = {RFF:.0%} x {FW_C:,.0f} = {W_res_C:,.0f} lb')
print(f'  W0 = MTOW                      = {W0_C:,} lb')
print(f'  W1 = OEW + Payload_C + Wres    = {OEW} + {payload_C:.0f} + {W_res_C:.0f} = {W1_C:,.0f} lb')
print(f'  Payload = {payload_C:,.0f} lb')
print(f'  Range   = {range_C:,.0f} nm')


# Point D, Ferry Range, Zero Payload
FW_D      = MFW
payload_D = 0
W_res_D   = RFF * FW_D
W0_D      = OEW + FW_D            # below MTOW -> no payload means lighter aircraft
W1_D      = OEW + W_res_D
range_D   = breguet(W0_D, W1_D)

print(f'\\nPoint D')
print(f'  Fuel loaded  = MFW             = {FW_D:,.0f} lb')
print(f'  W0 = OEW + MFW                 = {OEW} + {FW_D:.0f} = {W0_D:,.0f} lb  (< MTOW = {MTOW:,})')
print(f'  Reserve      = {RFF:.0%} x {FW_D:,.0f} = {W_res_D:,.0f} lb')
print(f'  W1 = OEW + Wres                = {OEW} + {W_res_D:.0f} = {W1_D:,.0f} lb')
print(f'  Payload = {payload_D} lb')
print(f'  Range   = {range_D:,.0f} nm')


# Point T, Design Mission
FW_T      = MTOW - OEW - DPW
W_res_T   = RFF * FW_T
payload_T = DPW
W0_T      = MTOW
W1_T      = OEW + DPW + W_res_T
range_T   = breguet(W0_T, W1_T)

print(f'\\nPoint T')
print(f'  Fuel loaded  = MTOW-OEW-DPW    = {MTOW}-{OEW}-{DPW} = {FW_T:,.1f} lb')
print(f'  Reserve      = {RFF:.0%} x {FW_T:,.1f} = {W_res_T:,.0f} lb')
print(f'  W0 = MTOW                       = {W0_T:,} lb')
print(f'  W1 = OEW + DPW + Wres           = {OEW} + {DPW} + {W_res_T:.0f} = {W1_T:,.0f} lb')
print(f'  Payload = {payload_T:,} lb')
print(f'  Range   = {range_T:,.0f} nm')


# Summary Table
summary = pd.DataFrame({
    'Point': ['A', 'B', 'C', 'D', 'T'],
    'Description': [
        'Zero range / max payload',
        'Harmonic range (max payload @ MTOW)',
        'Max-fuel range (full tanks @ MTOW)',
        'Ferry range (zero payload)',
        'Design / target mission',
    ],
    'Payload (lb)': [int(payload_A), int(payload_B), int(round(payload_C)),
                     int(payload_D), int(payload_T)],
    'TOW (lb)':     [int(OEW+payload_A), int(W0_B), int(W0_C),
                     int(round(W0_D)), int(W0_T)],
    'Fuel (lb)':    [0, int(FW_B), int(round(FW_C)), int(round(FW_D)), int(FW_T)],
    'W0 (lb)':      [int(OEW+payload_A), int(W0_B), int(W0_C),
                     int(round(W0_D)), int(W0_T)],
    'W1 (lb)':      [int(OEW+payload_A), int(round(W1_B)), int(round(W1_C)),
                     int(round(W1_D)), int(round(W1_T))],
    'Range (nm)':   [int(range_A), int(round(range_B)), int(round(range_C)),
                     int(round(range_D)), int(round(range_T))],
})
print("\\n", summary)


# Payload-Range Chart 
fig, ax = plt.subplots(figsize=(11, 7))


env_r = [range_A, range_B, range_C, range_D]
env_p = [payload_A, payload_B, payload_C, payload_D]

# Plot main envelope
ax.plot(env_r, env_p, color='#1f77b4', lw=3, marker='o', mfc='white', mew=2, label='Structural Envelope')

# Plot Design Point T
ax.scatter(range_T, payload_T, color='red', s=120, edgecolors='black', zorder=5, label=f'Design Point (T), 1836 lb @ 1927 nm')

# Clean Annotations (Offsets adjusted for no-overlap)
ann_data = [
    ('A', range_A, payload_A, (10, 10)),
    ('B', range_B, payload_B, (10, 10)),
    ('C', range_C, payload_C, (10, -20)),
    ('D', range_D, payload_D, (-60, 10)),
    ('T', range_T, payload_T, (-10, 15))
]

for label, x, y, offset in ann_data:
    ax.annotate(f'{label}: {int(x):,} nm', xy=(x, y), xytext=offset, 
                textcoords='offset points', fontweight='bold', fontsize=10)

# Aesthetics
ax.set_title('Payload-Range Chart', fontsize=16, pad=20)
ax.set_xlabel('Range (Nautical Miles)', fontsize=12)
ax.set_ylabel('Payload Weight (lb)', fontsize=12)
ax.grid(True, which='major', linestyle='--', alpha=0.5)
ax.set_ylim(0, MPW * 1.3)
ax.set_xlim(0, range_D * 1.1)
ax.legend(loc='upper right', frameon=True)

plt.tight_layout()
plt.show()