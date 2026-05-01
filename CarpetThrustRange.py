import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

# =============================================================
# Fighter Jet Carpet Plot - Altitude Labels + Clean Legend
# =============================================================

# ====================== USER INPUT SECTION ======================

rho_values = np.array([0.001756, 0.001496, 0.001267, 0.001066, 0.000891, 0.000738])

# Corresponding altitudes (10k to 35k ft, 5k increments)
altitudes = np.array([10000, 15000, 20000, 25000, 30000, 35000])   # feet

CD0_values = np.array([0.02304512, 0.0288064, 0.036008, 0.04501, 0.054012, 0.0648144])

# Fixed parameters
V_fps      = 850.0
Sref       = 600.0
W          = 52390.0
AR         = 2.028
e          = 0.82
k          = 1 / (np.pi * AR * e)
ct         = 0.00022222
CL         = 0.64180
CD_trim    = 0.0005
CD_wave    = 0.0005
CDi_const  = 0.0005
W_initial  = 52390.0
W_final    = 42540.0

# =============================================================

fig, ax = plt.subplots(figsize=(11, 8.5))   # Slightly taller for top legend

blue_cmap = get_cmap('Blues')
red_cmap = get_cmap('Reds')

# === Constant Altitude (rho) lines - Blue ===
for i, (rho, alt) in enumerate(zip(rho_values, altitudes)):
    x_list = []
    y_list = []
    color = blue_cmap(0.35 + 0.65 * i / (len(rho_values)-1))
    
    for CD0 in CD0_values:
        T_req = (0.5 * rho * V_fps**2 * CD0 * Sref) + \
                (2 * k * W**2) / (rho * V_fps**2 * Sref)
        
        CD_total = CD0 + CD_trim + CD_wave + CDi_const
        Range_ft = 2 * np.sqrt(2 / (rho * Sref)) * \
                   (1 / ct) * \
                   (np.sqrt(CL) / CD_total) * \
                   (np.sqrt(W_initial) - np.sqrt(W_final))
        
        Range_nm = Range_ft / 6076.12
        
        x_list.append(Range_nm)
        y_list.append(T_req)
    
    ax.plot(x_list, y_list, '-', linewidth=2.6, color=color,
            label=f'{alt:,} ft')

# === Constant CD0 lines - Red dashed ===
for i, CD0 in enumerate(CD0_values):
    x_list = []
    y_list = []
    color = red_cmap(0.35 + 0.65 * i / (len(CD0_values)-1))
    
    for rho in rho_values:
        T_req = (0.5 * rho * V_fps**2 * CD0 * Sref) + \
                (2 * k * W**2) / (rho * V_fps**2 * Sref)
        
        CD_total = CD0 + CD_trim + CD_wave + CDi_const
        Range_ft = 2 * np.sqrt(2 / (rho * Sref)) * \
                   (1 / ct) * \
                   (np.sqrt(CL) / CD_total) * \
                   (np.sqrt(W_initial) - np.sqrt(W_final))
        
        Range_nm = Range_ft / 6076.12
        
        x_list.append(Range_nm)
        y_list.append(T_req)
    
    ax.plot(x_list, y_list, '--', linewidth=2.2, color=color,
            label=f'CD₀ = {CD0:.5f}')

# ====================== Formatting ======================
ax.set_xlabel('Cruise Range [nautical miles]', fontsize=14)
ax.set_ylabel('Cruise Thrust Required [lbf]', fontsize=14)
ax.set_title('Carpet Plot\n'
             'CD0 and Cruise Altitude Effects on Required Cruise Thrust and Cruise Range', 
             fontsize=15, pad=25)

ax.grid(True, linestyle='--', alpha=0.6)

# Legend 
ax.legend(title='Constant Altitude (Blue)          Constant C$_{D0}$ (Red Dashed)',
          title_fontsize=11.5,
          fontsize=10.8, 
          loc='upper right', 
          ncol=2,
          frameon=True,
          edgecolor='gray')
plt.tight_layout()
plt.show()