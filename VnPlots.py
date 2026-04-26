import numpy as np
import matplotlib.pyplot as plt

def compute_stall_boundary(V_min, V_max, n_points, rho, CL_max, Weight, S_ref):
    print("Computing stall boundary...")
    stall_coeff = 0.5 * rho * CL_max / (Weight / S_ref)
    print(f"Wing loading W/S: {Weight / S_ref:.3f} lb/ft²")
    print(f"Stall coefficient: {stall_coeff:.8f}")
    V = np.linspace(V_min, V_max, n_points)
    n_Stall = stall_coeff * V**2
    return V, n_Stall, stall_coeff

# ================== Aircraft Parameters ==================
T_sl = 43000      # lbf (F135 sea-level static)
rho = 0.0008893   # slug/ft³ at 30,000 ft
rho_sl = 0.002377 # slug/ft³ sea level
T = T_sl * (rho / rho_sl)   # altitude-corrected thrust

W = 54748.05      # lbf
S = 600           # ft²
CD0 = 0.007       # parasite drag coefficient
e = 0.5
AR = 2.028
k = 1 / (np.pi * e * AR)

V_min = 0.0
n_points = 200

a_sound = 994.7   # ft/s at 30,000 ft

print(f"Thrust at 30,000 ft: {T:.1f} lbf")

# ================== Wave Drag Models ==================
Amax = 63.0
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

def drag(V):
    CL = 2 * W / (rho * S * V**2)
    CD = CD0 + k * CL**2 + CD_wave(mach(V))
    return 0.5 * rho * S * V**2 * CD

def f(V):
    return drag(V) - T

# ================== Find V_max ==================
V_grid = np.linspace(200, 5000, 20000)
idx = np.where(np.diff(np.sign(f(V_grid))) != 0)[0]
V_roots = []
for i in idx:
    v0, v1 = V_grid[i], V_grid[i+1]
    root = v0 - f(v0) * (v1 - v0) / (f(v1) - f(v0))
    V_roots.append(root)

V_max = max(V_roots)
print(f"\nMaximum level-flight speed: {V_max:.1f} ft/s  ({V_max/1.6878:.1f} knots,  Mach {V_max/a_sound:.2f})")

# ================== Stall Boundary (Fixed) ==================
CL_max = 1.55          # ← Realistic value for clean 5th-gen fighter
Boundary_V, Boundary_n_Stall, stall_coeff = compute_stall_boundary(
    V_min, V_max, n_points, rho, CL_max, W, S
)

n_pos_limit = 7.5
n_neg_limit = -4.5

V_corner_pos = np.sqrt(n_pos_limit / stall_coeff)
V_corner_neg = np.sqrt(abs(n_neg_limit) / stall_coeff)

print(f"\nPositive corner velocity (+{n_pos_limit} g): {V_corner_pos:.1f} ft/s ({V_corner_pos/1.6878:.1f} knots)")
print(f"Negative corner velocity ({n_neg_limit} g): {V_corner_neg:.1f} ft/s ({V_corner_neg/1.6878:.1f} knots)")

# ================== Plotting Data ==================
v_stall_pos = np.linspace(0, V_corner_pos, n_points)
n_stall_pos = stall_coeff * v_stall_pos**2

v_stall_neg = np.linspace(0, V_corner_neg, n_points)
n_stall_neg = -stall_coeff * v_stall_neg**2

v_limit_pos = np.linspace(V_corner_pos, V_max, n_points)
v_limit_neg = np.linspace(V_corner_neg, V_max, n_points)

v_max_line = np.full(n_points, V_max)
n_vertical = np.linspace(n_neg_limit, n_pos_limit, n_points)

# ================== Final V-n Diagram ==================
plt.figure(figsize=(16, 9))
plt.title('V-n Diagram: 5th Generation Fighter Maneuvering Envelope (30,000 ft)', fontsize=16)
plt.xlabel('True Airspeed (ft/s)', fontsize=12)
plt.ylabel('Load Factor n (g)', fontsize=12)

# Stall boundaries
plt.plot(v_stall_pos, n_stall_pos, 'b-', linewidth=3, label='Stall Boundary (Positive)')
plt.plot(v_stall_neg, n_stall_neg, 'b-', linewidth=3)

# Structural limits
plt.plot(v_limit_pos, np.full_like(v_limit_pos, n_pos_limit), 'r-', linewidth=3, label=f'Structural Limit (+{n_pos_limit} g)')
plt.plot(v_limit_neg, np.full_like(v_limit_neg, n_neg_limit), 'r-', linewidth=3, label=f'Structural Limit ({n_neg_limit} g)')

# Maximum speed limit
plt.plot(v_max_line, n_vertical, 'g-', linewidth=3, label=f'Max Speed Limit ({V_max/1.6878:.0f} knots)')

# Shaded safe envelope
plt.fill_between(v_stall_pos, 0, n_stall_pos, alpha=0.12, color='blue')
plt.fill_between(v_stall_neg, n_stall_neg, 0, alpha=0.12, color='blue')
plt.fill_between(v_limit_pos, n_stall_pos[-1], n_pos_limit, alpha=0.12, color='blue')
plt.fill_between(v_limit_neg, n_neg_limit, n_stall_neg[-1], alpha=0.12, color='blue')

plt.grid(True, alpha=0.4)
plt.legend(loc='best', fontsize=11)
plt.xlim(0, V_max * 1.02)
plt.ylim(n_neg_limit - 0.5, n_pos_limit + 0.5)
plt.show()