import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.optimize import root_scalar

# ================== Aircraft Parameters ==================
T_sl_max = 43000.0
delta_T = 1
K_T = 0.75

rho = 0.0008893
rho_sl = 0.002377
a_sound = 994.7

S = 600
CD0 = 0.007
e = 0.5
AR = 2.028
k = 1 / (np.pi * e * AR)

CL_max = 1.96
n_pos_limit = 8
n_neg_limit = -4

# ================== Wave Drag Models ==================
Amax = 54
length = 47.78
Ewd = 1.4
M_dd = 0.9

def mach(V):
    return V / a_sound

def sears_haack_drag():
    return (9 * np.pi / 2.0) * (Amax / length)**2

def wave_drag(M):
    M = np.asarray(M, dtype=float)
    D_q_SH = sears_haack_drag()
    result = np.zeros_like(M)
    mask = M >= 1.2
    if np.any(mask):
        mach_corr = 1 - 0.386 * (M[mask] - 1.2)**0.57
        result[mask] = Ewd * mach_corr * D_q_SH / S
    return result

def transonic_drag_rise(M):
    M = np.asarray(M, dtype=float)
    result = np.zeros_like(M)
    mid = (M >= M_dd) & (M < 1.2)
    if np.any(mid):
        t = (M[mid] - M_dd) / (1.2 - M_dd)
        result[mid] = 0.01 * t**3
    result[M >= 1.2] = 0.01
    return result

def CD_wave(M):
    return wave_drag(M) + transonic_drag_rise(M)

# ================== Thrust Model ==================
def get_thrust(V):
    M = mach(V)
    density_ratio = rho / rho_sl
    return T_sl_max * delta_T * density_ratio * (1 + K_T * M)

# ================== Drag Model ==================
def drag(V, Weight):
    CL = 2 * Weight / (rho * S * V**2)
    CD = CD0 + k * CL**2 + CD_wave(mach(V))
    q = 0.5 * rho * V**2
    return q * S * CD

def excess_thrust(V, Weight):
    return drag(V, Weight) - get_thrust(V)

# ================== Vmax Finder ==================
def find_vmax(Weight):
    def f(V):
        return excess_thrust(V, Weight)

    V_grid = np.linspace(800, 4000, 6000)
    f_grid = np.array([f(v) for v in V_grid])

    sign_changes = np.where(np.diff(np.sign(f_grid)) != 0)[0]

    if len(sign_changes) > 0:
        i = sign_changes[-1]
        bracket = [V_grid[max(0, i-20)], V_grid[min(len(V_grid)-1, i+30)]]
    else:
        bracket = [1100, 2500]

    try:
        sol = root_scalar(f, bracket=bracket, method='brentq', xtol=1e-6)
        if sol.converged:
            return sol.root
    except:
        pass

    idx = np.argmin(np.abs(f_grid))
    return V_grid[idx]

# ================== V-n Plot Function ==================
def plot_vn(Weight):
    print(f"\n=== Weight: {Weight:.0f} lb ===")

    V_C = find_vmax(Weight)
    V_D = 1.13 * V_C

    print(f"V_C (Vmax): {V_C:.1f} ft/s ({V_C/1.6878:.1f} knots, Mach {mach(V_C):.3f})")
    print(f"V_D (Dive): {V_D:.1f} ft/s ({V_D/1.6878:.1f} knots, Mach {mach(V_D):.3f})")

    k = 0.5 * rho * CL_max / (Weight / S)
    V_A_pos = np.sqrt(n_pos_limit / k)
    V_A_neg = np.sqrt(abs(n_neg_limit) / k)
    V_1g    = np.sqrt(1 / k)

    q_max = 0.5 * rho * V_D**2

    n = 400

    v_stall_pos = np.linspace(1, V_A_pos, n)
    n_stall_pos = k * v_stall_pos**2

    v_stall_neg = np.linspace(1, V_A_neg, n)
    n_stall_neg = -k * v_stall_neg**2

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_title(f'V-n Diagram at 30,000 ft  |  Weight = {Weight:,.0f} lb',
                 fontsize=14, fontweight='bold', pad=14)
    ax.set_xlabel('True Airspeed (ft/s)', fontsize=12)
    ax.set_ylabel('Load Factor  n  (g)', fontsize=12)

    ax.plot(v_stall_pos, n_stall_pos, color='royalblue', lw=2.5, label='Positive Stall')
    ax.plot(v_stall_neg, n_stall_neg, color='deepskyblue', lw=2.5, label='Negative Stall')

    ax.plot([V_A_pos, V_D], [n_pos_limit, n_pos_limit],
            color='crimson', lw=2.5, label=f'Positive Load Factor Limit ({n_pos_limit} g)')

    ax.plot([V_A_neg, V_C], [n_neg_limit, n_neg_limit],
            color='firebrick', lw=2.5, label=f'Negative Load Factor Limit ({n_neg_limit} g)')

    ax.plot([V_C, V_C], [n_neg_limit, n_pos_limit],
            color='seagreen', lw=2.5,
            label=f'Max Level Flight Speed ($V_C$) = {V_C:.0f} ft/s')

    ax.plot([V_D, V_D], [0, n_pos_limit],
            color='mediumpurple', lw=2.5,
            label=f'Dive Speed ($V_D$) = {V_D:.0f} ft/s')

    ax.plot([V_C, V_D], [n_neg_limit, 0], color='k', lw=2.5, label='Dive Boundary')

    ax.axhline(1, ls='--', color='gray', lw=1, alpha=0.7)
    ax.axhline(-1, ls='--', color='gray', lw=1, alpha=0.7)

    ax.fill_between(v_stall_pos, 0, n_stall_pos, alpha=0.07, color='royalblue')
    ax.fill_between([V_A_pos, V_D], n_stall_pos[-1], n_pos_limit, alpha=0.07, color='royalblue')
    ax.fill_between(v_stall_neg, n_stall_neg, 0, alpha=0.07, color='deepskyblue')
    ax.fill_between([V_A_neg, V_C], n_neg_limit, n_stall_neg[-1], alpha=0.07, color='deepskyblue')

    ax.scatter(V_A_pos, n_pos_limit, color='royalblue', s=70, zorder=5)
    ax.scatter(V_A_neg, n_neg_limit, color='firebrick', s=70, zorder=5)
    ax.scatter(V_1g, 1, color='gray', s=50, zorder=5)
    ax.scatter(V_1g, -1, color='gray', s=50, zorder=5)
    ax.scatter(V_D, n_pos_limit, color='black', s=90, zorder=5)

    def label(x, y, txt, ha='center', va='center', color='black', fontsize=9.5, offset=(0, 0)):
        ax.text(x + offset[0], y + offset[1], txt,
                ha=ha, va=va, fontsize=fontsize, color=color,
                path_effects=[pe.withStroke(linewidth=3, foreground='white')])

    label(V_A_pos, n_pos_limit,
          f'Positive Cornering Speed ({V_A_pos:.0f} ft/s)',
          ha='center', va='bottom', color='royalblue', offset=(0, 0.15))

    label(V_A_neg, n_neg_limit,
          f'Negative Cornering Speed ({V_A_neg:.0f} ft/s)',
          ha='center', va='top', color='firebrick', offset=(0, -0.15))

    label(V_1g, 1,
          f'$V_{{+1g}}$ = {V_1g:.0f} ft/s',
          ha='left', va='bottom', color='dimgray', offset=(12, 0.12))

    label(V_1g, -1,
          f'$V_{{-1g}}$ = {V_1g:.0f} ft/s',
          ha='left', va='top', color='dimgray', offset=(12, -0.12))

    label(V_D, n_pos_limit,
          f'$q_{{max}}$ = {q_max:.0f} lb/ft²',
          ha='left', va='bottom', color='black', offset=(8, 0.15))

    x_right = V_D * 1.03
    ax.text(x_right, 1.18, '+1g', ha='left', va='center', fontsize=9, color='gray')
    ax.text(x_right, -1.18, '−1g', ha='left', va='center', fontsize=9, color='gray')

    ax.grid(True, alpha=0.35)
    ax.legend(loc='upper left', fontsize=9.5, framealpha=0.9)
    ax.set_xlim(0, V_D * 1.07)
    ax.set_ylim(n_neg_limit - 1.2, n_pos_limit + 1.0)
    fig.tight_layout()
    plt.show()

# ================== Run Both Cases ==================
weights = [54748.05, 31862.19]

for W in weights:
    plot_vn(W)

# ================== Weight Sensitivity Test ==================
print("\n" + "="*70)
print("WEIGHT SENSITIVITY TEST")
print("="*70)

test_weights = [31862.19, 54748.05, 65000]
for test_W in test_weights:
    test_Vmax = find_vmax(test_W)
    if test_Vmax:
        test_k = 0.5 * rho * CL_max / (test_W / S)
        test_corner = np.sqrt(n_pos_limit / test_k)
        print(f"W = {test_W:6.0f} lb → V_max = {test_Vmax:.1f} ft/s ({test_Vmax/1.6878:.1f} knots) → Corner V = {test_corner:.1f} ft/s")