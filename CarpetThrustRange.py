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

CD0_values = np.array([0.00645165, 0.0071685, 0.007965, 0.00885, 0.009735, 0.0107085])
# Fixed parameters
V_fps      = 800
Sref       = 600.0
W          = 52390.0
AR         = 2.028
e          = 0.82
k          = 1 / (np.pi * AR * e)
ct         = 0.00022222
CL         = 0.64180
CD         = 0.0386
CD0Design  = 0.00885
W_initial  = 52390.0
W_final    = 40000.0

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
        
        CD_total = (CD0 - CD0Design) + CD
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
        
        CD_total = (CD0 - CD0Design) + CD
        Range_ft = 2 * np.sqrt(2 / (rho * Sref)) * \
                   (1 / ct) * \
                   (np.sqrt(CL) / CD_total) * \
                   (np.sqrt(W_initial) - np.sqrt(W_final))
        
        Range_nm = Range_ft / 6076.12
        
        x_list.append(Range_nm)
        y_list.append(T_req)
    
    ax.plot(x_list, y_list, '--', linewidth=2.2, color=color,
            label=f'CD₀ = {CD0:.5f}')

# ====================== PLOTTING SECTION ======================

fig, ax = plt.subplots(figsize=(12.5, 9))   

blue_cmap = get_cmap('Blues')
red_cmap = get_cmap('Reds')

all_x = []
all_y = []

# === Constant Altitude (rho) lines - Blue ===
for i, (rho, alt) in enumerate(zip(rho_values, altitudes)):
    x_list = []
    y_list = []
    color = blue_cmap(0.35 + 0.65 * i / (len(rho_values)-1))
    
    for CD0 in CD0_values:
        T_req = (0.5 * rho * V_fps**2 * CD0 * Sref) + \
                (2 * k * W**2) / (rho * V_fps**2 * Sref)
        
        CD_total = (CD0 - CD0Design) + CD
        Range_ft = 2 * np.sqrt(2 / (rho * Sref)) * \
                   (1 / ct) * \
                   (np.sqrt(CL) / CD_total) * \
                   (np.sqrt(W_initial) - np.sqrt(W_final))
        
        Range_nm = Range_ft / 6076.12
        
        x_list.append(Range_nm)
        y_list.append(T_req)
        all_x.append(Range_nm)
        all_y.append(T_req)
    
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
        
        CD_total = (CD0 - CD0Design) + CD
        Range_ft = 2 * np.sqrt(2 / (rho * Sref)) * \
                   (1 / ct) * \
                   (np.sqrt(CL) / CD_total) * \
                   (np.sqrt(W_initial) - np.sqrt(W_final))
        
        Range_nm = Range_ft / 6076.12
        
        x_list.append(Range_nm)
        y_list.append(T_req)
        all_x.append(Range_nm)
        all_y.append(T_req)
    
    ax.plot(x_list, y_list, '--', linewidth=2.2, color=color,
            label=f'CD₀ = {CD0:.5f}')

# ====================== PLOTTING SECTION ======================

fig, ax = plt.subplots(figsize=(12.5, 9))   

blue_cmap = get_cmap('Blues')
red_cmap = get_cmap('Reds')

all_x = []
all_y = []

# === Constant Altitude (rho) lines - Blue ===
for i, (rho, alt) in enumerate(zip(rho_values, altitudes)):
    x_list = []
    y_list = []
    color = blue_cmap(0.35 + 0.65 * i / (len(rho_values)-1))
    
    for CD0 in CD0_values:
        T_req = (0.5 * rho * V_fps**2 * CD0 * Sref) + \
                (2 * k * W**2) / (rho * V_fps**2 * Sref)
        
        CD_total = (CD0 - CD0Design) + CD
        Range_ft = 2 * np.sqrt(2 / (rho * Sref)) * \
                   (1 / ct) * \
                   (np.sqrt(CL) / CD_total) * \
                   (np.sqrt(W_initial) - np.sqrt(W_final))
        
        Range_nm = Range_ft / 6076.12
        
        x_list.append(Range_nm)
        y_list.append(T_req)
        all_x.append(Range_nm)
        all_y.append(T_req)
    
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
        
        CD_total = (CD0 - CD0Design) + CD
        Range_ft = 2 * np.sqrt(2 / (rho * Sref)) * \
                   (1 / ct) * \
                   (np.sqrt(CL) / CD_total) * \
                   (np.sqrt(W_initial) - np.sqrt(W_final))
        
        Range_nm = Range_ft / 6076.12
        
        x_list.append(Range_nm)
        y_list.append(T_req)
        all_x.append(Range_nm)
        all_y.append(T_req)
    
    ax.plot(x_list, y_list, '--', linewidth=2.2, color=color,
            label=f'CD₀ = {CD0:.5f}')

# ====================== DESIGN POINT (30,000 ft) ======================
design_rho = rho_values[4]        # ← Correct index for 30,000 ft
design_cd0 = 0.00885

CD_total_design = (design_cd0 - CD0Design) + CD

T_req_design = (0.5 * design_rho * V_fps**2 * design_cd0 * Sref) + \
               (2 * k * W**2) / (design_rho * V_fps**2 * Sref)

Range_ft_design = 2 * np.sqrt(2 / (design_rho * Sref)) * \
                  (1 / ct) * \
                  (np.sqrt(CL) / CD_total_design) * \
                  (np.sqrt(W_initial) - np.sqrt(W_final))

Range_nm_design = Range_ft_design / 6076.12

# Plot design point
ax.plot(Range_nm_design, T_req_design, 'o', color='gold', markersize=14, 
        markeredgecolor='black', markeredgewidth=2.5, zorder=10, 
        label='_nolegend_')

ax.annotate('Current Design Point\n(30,000 ft • CD₀ = 0.00885)', 
            xy=(Range_nm_design, T_req_design),
            xytext=(Range_nm_design -200, T_req_design + 850),
            fontsize=11.5,
            fontweight='bold',
            ha='right',
            arrowprops=dict(arrowstyle='->', color='black', lw=1.8),
            bbox=dict(boxstyle="round,pad=0.5", facecolor="yellow", alpha=0.95))

# ====================== Formatting ======================
ax.set_xlabel('Cruise Range [nautical miles]', fontsize=14)
ax.set_ylabel('Cruise Thrust Required [lbf]', fontsize=14)
ax.set_title('Carpet Plot\n'
             'CD0 and Cruise Altitude Effects on Required Cruise Thrust and Cruise Range', 
             fontsize=15, pad=25)

ax.grid(True, linestyle='--', alpha=0.6)

# Gentle zoom out
ax.set_xlim(left=min(all_x)*0.93, right=max(all_x)*1.08)
ax.set_ylim(bottom=min(all_y)*0.88, top=max(all_y)*1.14)

# Legend 
ax.legend(title='Constant Altitude (Blue)          Constant C$_{D0}$ (Red Dashed)',
          title_fontsize=11.5,
          fontsize=10.5, 
          loc='upper right', 
          ncol=2,
          frameon=True,
          edgecolor='gray')

plt.tight_layout()
plt.show()