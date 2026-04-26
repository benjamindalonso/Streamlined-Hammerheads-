import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

def compute_stall_boundary(V_min, V_max, n_points, rho, CL_max, Weight, S_ref):
    print("Computing stall boundary...")
    stall_coeff = 0.5 * rho * CL_max / (Weight / S_ref)
    print(f"W/S: {Weight / S_ref:.3f} N/m^2")
    print(f"Stall coefficient: {stall_coeff:.6f}")

    V = np.linspace(V_min, V_max, n_points)
    n_Stall = stall_coeff * V**2

    return V, n_Stall, stall_coeff

# Given parameters
T_sl = 43000      # lbf, F135 sea-level static thrust
rho = 0.0008893    # slug/ft^3 at 30,000 ft
rho_sl = 0.002377 # slug/ft^3 at sea level
T = T_sl * (rho / rho_sl)
W =  54748.05     # lbf
S = 600           # ft^2
CD0 = 0.007       # from OpenVSP parasite drag (clean)
n_points = 100     # number of points in the speed range
e = 0.5
AR = 2.028
k = 1 / (np.pi * e * AR)
V_min = 0          # ft/s, minimum speed

print(f"Using altitude-corrected thrust T = {T:.1f} lbf at 30,000 ft")

# Speed of sound at 30,000 ft (ft/s)
a_sound = 994.7

# Additional wave-drag model constants from A2_MachEstimatePlot
Amax = 63.0        # ft^2 reference max cross-sectional area
length = 47.78     # ft
Ewd = 1.4          # wave drag efficiency factor
M_dd = 0.9         # drag divergence Mach


def mach(V, a=a_sound):
    return V / a


def sears_haack_drag(Amax, length):
    return (9 * np.pi / 2.0) * (Amax / length)**2


def wave_drag(M):
    M = np.asarray(M, dtype=float)
    D_q_SH = sears_haack_drag(Amax, length)
    
    # Initialize output
    result = np.zeros_like(M)
    
    # Only compute for M >= 1.2
    mask = M >= 1.2
    mach_corr = 1 - 0.386 * (M[mask] - 1.2)**0.57
    result[mask] = Ewd * mach_corr * D_q_SH / S
    
    return float(result) if result.shape == () else result


def transonic_drag_rise(M):
    M = np.asarray(M, dtype=float)
    result = np.zeros_like(M)
    
    # M < M_dd: return 0.0 (already initialized to 0)
    
    # M_dd <= M < 1.2: cubic rise
    mid_mask = (M >= M_dd) & (M < 1.2)
    t = (M[mid_mask] - M_dd) / (1.2 - M_dd)
    result[mid_mask] = 0.01 * t**3
    
    # M >= 1.2: return 0.01
    high_mask = M >= 1.2
    result[high_mask] = 0.01
    
    return float(result) if result.shape == () else result


def CD_wave(M):
    M = np.asarray(M, dtype=float)
    return wave_drag(M) + transonic_drag_rise(M)


def drag(V):
    CL = 2 * W / (rho * S * V**2)
    CD = CD0 + k * CL**2 + CD_wave(mach(V))
    return 0.5 * rho * S * V**2 * CD


def f(V):
    return drag(V) - T

# Find the high-speed root using a simple bracket search
V_min_HighSpeed = 200.0
V_max_HighSpeed = 5000.0
V_grid = np.linspace(V_min_HighSpeed, V_max_HighSpeed, 10000)
F = f(V_grid)

# Find sign changes
idx = np.where(np.sign(F[:-1]) != np.sign(F[1:]))[0]
if len(idx) == 0:
    raise RuntimeError("No root found for drag = thrust with wave drag included.")

# Choose the highest-speed root
root_intervals = [(V_grid[i], V_grid[i+1]) for i in idx]
V_roots = []
for v0, v1 in root_intervals:
    f0 = f(v0)
    f1 = f(v1)
    V_roots.append(v0 - f0 * (v1 - v0) / (f1 - f0))

V_max = max(V_roots)
print(f"Upper root V = {V_max:.3f} ft/s")

Vtest = V_max
CL = 2*W/(rho*S*Vtest**2)
Di = 0.5*rho*S*Vtest**2 * k * CL**2
D0 = 0.5*rho*S*Vtest**2 * CD0
Dw = 0.5*rho*S*Vtest**2 * CD_wave(mach(Vtest))

print("Induced:", Di)
print("Parasite:", D0)
print("Wave:", Dw)

# ============================================================================
# Calculate Stall Boundary and Initialize Velocity Array
# ============================================================================
Boundary_V, Boundary_n_Stall, stall_coeff = compute_stall_boundary(
    V_min, V_max, n_points, rho, CL, W, S
)
velocity_range = np.linspace(V_min, V_max, n_points)

# ============================================================================
# Plot 1: Basic Stall Boundary
# ============================================================================
plt.figure(figsize=(16, 9))
plt.title('V-n Diagram: Stall Boundary Only')
plt.xlabel('Velocity (ft/s)')
plt.ylabel('Load Factor n (-)')
plt.plot(Boundary_V, Boundary_n_Stall, 'b-', linewidth=2, label='Stall (Positive)')
plt.plot(Boundary_V, -Boundary_n_Stall, 'b-', linewidth=2, label='Stall (Negative)')
plt.grid(True, alpha=0.3)
plt.legend(loc='best')
plt.show()

# ============================================================================
# Define Load Factors for 5th Generation Fighter
# ============================================================================
n_pos_limit = np.ones(n_points) * 7.5   # Positive limit load factor
n_neg_limit = np.ones(n_points) * -4.5  # Negative limit load factor
print(f"Positive limit load factor: {n_pos_limit[0]:.2f} g")
print(f"Negative limit load factor: {n_neg_limit[0]:.2f} g")

# ============================================================================
# Plot 2: Stall Boundary with Limit Load Lines
# ============================================================================
plt.figure(figsize=(16, 9))
plt.title('V-n Diagram: Stall Boundary + Limit Load Factors')
plt.xlabel('Velocity (ft/s)')
plt.ylabel('Load Factor n (-)')
plt.plot(Boundary_V, Boundary_n_Stall, 'b-', linewidth=2, label='Stall (Positive)')
plt.plot(Boundary_V, -Boundary_n_Stall, 'b-', linewidth=2, label='Stall (Negative)')
plt.plot(velocity_range, n_pos_limit, 'r--', linewidth=2, label=f'Limit Load: +{n_pos_limit[0]:.1f} g')
plt.plot(velocity_range, n_neg_limit, 'r--', linewidth=2, label=f'Limit Load: {n_neg_limit[0]:.1f} g')
plt.grid(True, alpha=0.3)
plt.legend(loc='best')
plt.show()

# ============================================================================
# Calculate Stall Envelope Corners (Intersection Points)
# ============================================================================
def find_corner_velocity(stall_coeff, load_factor):
    """Find velocity where stall boundary intersects a load factor limit."""
    corner_speed = np.sqrt(abs(load_factor) / stall_coeff)
    print(f"Corner velocity at n = {load_factor:+.2f}: {corner_speed:.2f} ft/s")
    return corner_speed

V_corner_pos = find_corner_velocity(stall_coeff, n_pos_limit[0])
V_corner_neg = find_corner_velocity(stall_coeff, n_neg_limit[0])

# Define stall boundary segments up to corner velocities
v_stall_pos_segment = np.linspace(0, V_corner_pos, n_points)
v_stall_neg_segment = np.linspace(0, V_corner_neg, n_points)

n_stall_pos_segment = stall_coeff * v_stall_pos_segment**2
n_stall_neg_segment = -stall_coeff * v_stall_neg_segment**2

# ============================================================================
# Plot 3: Stall Boundary Segments with Corners
# ============================================================================
plt.figure(figsize=(16, 9))
plt.title('V-n Diagram: Stall Boundary Segments at Corner Velocities')
plt.xlabel('Velocity (ft/s)')
plt.ylabel('Load Factor n (-)')
plt.plot(v_stall_pos_segment, n_stall_pos_segment, 'b-', linewidth=2.5, label='Stall Boundary (Positive)')
plt.plot(v_stall_neg_segment, n_stall_neg_segment, 'b-', linewidth=2.5, label='Stall Boundary (Negative)')
plt.scatter([V_corner_pos], [n_pos_limit[0]], color='red', s=100, zorder=5, label='Corner Point (Positive)')
plt.scatter([V_corner_neg], [n_neg_limit[0]], color='red', s=100, zorder=5, label='Corner Point (Negative)')
plt.grid(True, alpha=0.3)
plt.legend(loc='best')
plt.show()

# ============================================================================
# Define Extended Limit Load Line Segments
# ============================================================================
v_limit_pos_extended = np.linspace(V_corner_pos, V_max, n_points)
v_limit_neg_extended = np.linspace(V_corner_neg, V_max, n_points)
n_limit_pos_extended = np.ones_like(v_limit_pos_extended) * n_pos_limit[0]
n_limit_neg_extended = np.ones_like(v_limit_neg_extended) * n_neg_limit[0]

# Define speed limit line at V_max
v_exceed_limit = np.ones(n_points) * V_max
n_exceed_range = np.linspace(n_neg_limit[0], n_pos_limit[0], n_points)

# ============================================================================
# Plot 4: Complete V-n Envelope (Stall Boundary + Extended Limits)
# ============================================================================
plt.figure(figsize=(16, 9))
plt.title('V-n Diagram: Complete Flight Envelope')
plt.xlabel('Velocity (ft/s)')
plt.ylabel('Load Factor n (-)')

# Stall boundaries
plt.plot(v_stall_pos_segment, n_stall_pos_segment, 'b-', linewidth=2.5, label='Stall Boundary')
plt.plot(v_stall_neg_segment, n_stall_neg_segment, 'b-', linewidth=2.5)

# Extended limit loads
plt.plot(v_limit_pos_extended, n_limit_pos_extended, 'r-', linewidth=2.5, label='Positive Limit Load')
plt.plot(v_limit_neg_extended, n_limit_neg_extended, 'r-', linewidth=2.5, label='Negative Limit Load')

# Speed limit
plt.plot(v_exceed_limit, n_exceed_range, 'g-', linewidth=2.5, label='Maximum Speed Limit')

plt.grid(True, alpha=0.3)
plt.legend(loc='best')
plt.show()

# ============================================================================
# Plot 5: Final Maneuvering Envelope with All Boundaries
# ============================================================================
plt.figure(figsize=(16, 9))
plt.title('V-n Diagram: Final Maneuvering Envelope for 5th Generation Fighter')
plt.xlabel('Velocity (ft/s)')
plt.ylabel('Load Factor n (-)')

# Stall boundaries (blue)
plt.plot(v_stall_pos_segment, n_stall_pos_segment, 'b-', linewidth=2.5, label='Stall Boundary')
plt.plot(v_stall_neg_segment, n_stall_neg_segment, 'b-', linewidth=2.5)

# Extended limit loads (red)
plt.plot(v_limit_pos_extended, n_limit_pos_extended, 'r-', linewidth=2.5, label='Structural Limit (Positive)')
plt.plot(v_limit_neg_extended, n_limit_neg_extended, 'orange', linewidth=2.5, label='Structural Limit (Negative)')

# Speed limit (green)
plt.plot(v_exceed_limit, n_exceed_range, 'g-', linewidth=2.5, label='Speed Limit')

# Shade the safe flight envelope
plt.fill_between(v_stall_pos_segment, 0, n_stall_pos_segment, alpha=0.1, color='blue')
plt.fill_between(v_stall_neg_segment, n_stall_neg_segment, 0, alpha=0.1, color='blue')
plt.fill_between(v_limit_pos_extended, n_stall_pos_segment[-1], n_limit_pos_extended, alpha=0.1, color='blue')
plt.fill_between(v_limit_neg_extended, n_limit_neg_extended, n_stall_neg_segment[-1], alpha=0.1, color='blue')

plt.grid(True, alpha=0.3)
plt.legend(loc='best', fontsize=10)
plt.xlim(0, V_max)
plt.ylim(min(n_neg_limit[0], -5), max(n_pos_limit[0], 8))
plt.show()
