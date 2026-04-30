import numpy as np
import matplotlib.pyplot as plt
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
n_pos_limit = 7.5
n_neg_limit = -4.5

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

    # fallback
    idx = np.argmin(np.abs(f_grid))
    return V_grid[idx]

# ================== V-n Plot Function ==================
def plot_vn(Weight):
    print(f"\n=== Weight: {Weight:.0f} lb ===")

    V_max = find_vmax(Weight)
    print(f"V_max: {V_max:.1f} ft/s ({V_max/1.6878:.1f} knots, Mach {mach(V_max):.3f})")

    stall_coeff = 0.5 * rho * CL_max / (Weight / S)
    V_corner_pos = np.sqrt(n_pos_limit / stall_coeff)
    V_corner_neg = np.sqrt(abs(n_neg_limit) / stall_coeff)

    n_points = 200

    v_stall_pos = np.linspace(0, V_corner_pos, n_points)
    n_stall_pos = stall_coeff * v_stall_pos**2

    v_stall_neg = np.linspace(0, V_corner_neg, n_points)
    n_stall_neg = -stall_coeff * v_stall_neg**2

    v_limit_pos = np.linspace(V_corner_pos, V_max, n_points)
    v_limit_neg = np.linspace(V_corner_neg, V_max, n_points)

    v_max_line = np.full(n_points, V_max)
    n_vertical = np.linspace(n_neg_limit, n_pos_limit, n_points)

    # ===== Plot =====
    plt.figure(figsize=(14, 8))
    plt.title(f'V-n Diagram at 30,000 ft (Weight = {Weight:.0f} lb)')
    plt.xlabel('True Airspeed (ft/s)')
    plt.ylabel('Load Factor n (g)')

    plt.plot(v_stall_pos, n_stall_pos, 'b-', linewidth=3, label='Stall (+)')
    plt.plot(v_stall_neg, n_stall_neg, 'b-', linewidth=3, label='Stall (-)')

    plt.plot(v_limit_pos, np.full_like(v_limit_pos, n_pos_limit), 'r-', linewidth=3,
             label=f'+{n_pos_limit} g limit')
    plt.plot(v_limit_neg, np.full_like(v_limit_neg, n_neg_limit), 'r-', linewidth=3,
             label=f'{n_neg_limit} g limit')

    plt.plot(v_max_line, n_vertical, 'g-', linewidth=3,
             label=f'Vmax ({V_max/1.6878:.0f} knots)')

    # Shaded region
    plt.fill_between(v_stall_pos, 0, n_stall_pos, alpha=0.1)
    plt.fill_between(v_limit_pos, n_stall_pos[-1], n_pos_limit, alpha=0.1)
    plt.fill_between(v_stall_neg, n_stall_neg, 0, alpha=0.1)
    plt.fill_between(v_limit_neg, n_neg_limit, n_stall_neg[-1], alpha=0.1)

    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.xlim(0, V_max * 1.05)
    plt.ylim(n_neg_limit - 0.5, n_pos_limit + 0.5)
    plt.tight_layout()
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
        test_stall_coeff = 0.5 * rho * CL_max / (test_W / S)
        test_corner = np.sqrt(n_pos_limit / test_stall_coeff)
        print(f"W = {test_W:6.0f} lb → V_max = {test_Vmax:.1f} ft/s ({test_Vmax/1.6878:.1f} knots) → Corner V = {test_corner:.1f} ft/s")