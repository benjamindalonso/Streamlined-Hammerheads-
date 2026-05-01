import numpy as np
import matplotlib.pyplot as plt

# ====================== AIRCRAFT INPUTS ======================
CD0 = 0.025      # ← Your known parasite drag coefficient
e = 0.85         # Oswald efficiency factor (tune this 0.75-0.9)
gamma = 1.4
p_sl = 2116.2    # Sea level pressure (lb/ft²). Lower for altitude.

# Ranges matching slide style
AR_values = np.array([6, 8, 10, 12, 14])           # Key ARs for constant Mach lines
Mach_values = np.array([0.60, 0.62, 0.64, 0.66, 0.68, 0.70])  # Key Mach for constant AR lines

num_points = 200   # For smooth curves
# ============================================================

def carpet_data(AR, Mach, CD0, e, gamma, p):
    """Compute (L/D)max and W/S"""
    CL_opt = np.sqrt(np.pi * AR * e * CD0)
    LD_max = 0.5 * np.sqrt(np.pi * AR * e / CD0)
    q = 0.5 * gamma * p * Mach**2
    WS = CL_opt * q
    return LD_max, WS

# =============== Generate Data ===============
fig, ax = plt.subplots(figsize=(11, 8))

# 1. Constant AR lines (vary Mach) - these are the curved lines in slides
for AR in AR_values:
    LDs = []
    WSs = []
    Mach_range = np.linspace(0.55, 0.75, num_points)
    for M in Mach_range:
        LD, WS = carpet_data(AR, M, CD0, e, gamma, p_sl)
        LDs.append(LD)
        WSs.append(WS)
    ax.plot(WSs, LDs, linewidth=2.2, label=f'AR = {AR}')

# 2. Constant Mach lines (vary AR) - diagonal-ish lines
for M in Mach_values:
    LDs = []
    WSs = []
    AR_range = np.linspace(5, 15, num_points)
    for AR in AR_range:
        LD, WS = carpet_data(AR, M, CD0, e, gamma, p_sl)
        LDs.append(LD)
        WSs.append(WS)
    ax.plot(WSs, LDs, '--', linewidth=1.8, alpha=0.85, label=f'M = {M:.2f}')

# ====================== Plot Formatting ======================
ax.set_xlabel('W/S (lb/ft²)', fontsize=12)
ax.set_ylabel('(L/D)$_{max}$', fontsize=12)
ax.set_title(f'Carpet Plot — (L/D)$_{{max}}$ vs Wing Loading\n'
             f'CD$_0$ = {CD0}, e = {e}, Sea Level', fontsize=14, pad=20)

ax.grid(True, alpha=0.3)
ax.set_xlim(50, 250)
ax.set_ylim(14, 22)

# Legend
ax.legend(loc='upper right', fontsize=10, ncol=2)

# Label the lines (like in the slides)
for AR in AR_values:
    # Label at high Mach end
    M_high = 0.72
    LD, WS = carpet_data(AR, M_high, CD0, e, gamma, p_sl)
    ax.text(WS*1.02, LD, f'AR={AR}', fontsize=10, va='center', fontweight='bold')

for M in Mach_values:
    AR_high = 14
    LD, WS = carpet_data(AR_high, M, CD0, e, gamma, p_sl)
    ax.text(WS*1.02, LD, f'M={M:.2f}', fontsize=10, va='center')

plt.tight_layout()
plt.show()

# Optional: Save data
import pandas as pd
results = []
for AR in AR_values:
    for M in Mach_values:
        LD, WS = carpet_data(AR, M, CD0, e, gamma, p_sl)
        results.append({'AR': AR, 'Mach': M, 'LD_max': round(LD,2), 'W/S': round(WS,1)})

