import numpy as np
import matplotlib.pyplot as plt

# Inputs
Sref = 600.0          # ft^2 
Amax = 63.0           # ft^2 (Ref. max cross-sectional area)
length = 47.78        # ft (aircraft length)
AR = 2.028            # aspect ratio
e = 0.8               # Oswald efficiency factor

CD0_subsonic = 0.007  # from OpenVSP parasite drag (clean)
CD_trim = 0.0         # from AVL trim drag estimate

Ewd = 1.4             # wave drag efficiency factor (1.4–2.0)
M_dd = 0.9           # drag divergence Mach

# for skin-friction scaling
Re_unit = 4.4e7       # Typical unit Reynolds number (per ft) for fighter at cruise
Re_L    = Re_unit * length   # Reynolds number based on aircraft length


def induced_drag(CL, AR, e):
    return CL**2 / (np.pi * AR * e)

def sears_haack_drag(Amax, length):
    return (9 * np.pi / 2.0) * (Amax / length)**2

def wave_drag(M):
    if M < 1.2:
        return 0.0
    
    D_q_SH = sears_haack_drag(Amax, length)
    
    # Mach correction (simplified from slides)
    mach_corr = (1 - 0.386 * (M - 1.2)**0.57)
    
    D_q_wave = Ewd * mach_corr * D_q_SH
    
    return D_q_wave / Sref

def transonic_drag_rise(M):
    if M < M_dd:
        return 0.0
    elif M < 1.2:
        # Cubic rise: smoother and more realistic than linear for plot
        t = (M - M_dd) / (1.2 - M_dd)
        return 0.01 * t**3          # reaches exactly 0.01 at M=1.2
    else:
        return 0.01

def calculate_cf(M, Re_L):
    """Raymer turbulent flat-plate skin-friction coefficient
    (Exact formula from the "Flat-Plate Skin-Friction Coefficient" slide)"""
    if Re_L <= 0:
        return 0.0
    term1 = (np.log10(Re_L))**2.58
    term2 = (1 + 0.144 * M**2)**0.65
    return 0.455 / (term1 * term2)

# Main Computation
Mach = np.linspace(0.2, 2.0, 100)
CD_total = []

# Reference Cf at low Mach (used to scale OpenVSP CD0_subsonic)
cf_ref = calculate_cf(0.1, Re_L)

# Assumed a representative lift coefficient for cruise (can be adjusted as needed)
CL = 0.0

for M in Mach:
    
    # Base parasite drag w/ skin friction scaling
    cf_current = calculate_cf(M, Re_L)
    CD0 = CD0_subsonic * (cf_current / cf_ref)
    
    # Induced drag
    CDi = induced_drag(CL, AR, e)
    
    # Transonic drag rise
    CD_trans = transonic_drag_rise(M)
    
    # Wave drag
    CD_wave = wave_drag(M)
    
    # Total drag 
    CD = CD0 + CDi + CD_trim + CD_trans + CD_wave
    
    CD_total.append(CD)

CD_total = np.array(CD_total)

# -----------------------------
# PLOTTING
# -----------------------------

plt.figure()
plt.plot(Mach, CD_total, label='Total Drag (Clean)')
plt.axvline(x=0.9, linestyle='--', label='Drag Divergence')
plt.axvline(x=1.2, linestyle='--', label='Full Supersonic Wave Drag')

# === Four RFP Points with annotations ===
rfp_points = [
    (0.85, "Strike SL Dash (req.)", 'green'),
    (0.90, "Strike SL Dash (des.)", 'darkgreen'),
    (1.60, "A/A Dash (req.)", 'darkorange'),
    (2.00, "A/A Dash (des.)", 'red')
]

for m, label, color in rfp_points:
    idx = np.argmin(np.abs(Mach - m))
    cd_val = CD_total[idx]

    # Only adjust the bottom (dark green) label
    if color == 'darkgreen':
        text_offset = (0, -30)   # move further down
        ha = 'center'
    else:
        text_offset = (8, 12)    # keep original
        ha = 'left'

    # Plot marker
    plt.plot(m, cd_val, marker='*', markersize=10, color=color)

    # Annotate
    plt.annotate(f'{cd_val:.4f}',
                 xy=(m, cd_val),
                 xytext=text_offset,
                 textcoords='offset pixels',
                 ha=ha,
                 fontsize=8,
                 color=color,
                 weight='bold')

plt.xlabel('Mach Number')
plt.ylabel('Drag Coefficient (CD)')
plt.title('High Mach Drag Estimate (Clean Configuration)')
plt.legend()
plt.grid()

plt.show()