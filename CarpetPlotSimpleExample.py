import numpy as np
import matplotlib.pyplot as plt

# =============================================================
# Simple carpet plot example
# =============================================================

# Two independent inputs (the "carpet parameters")
A_values = np.array([1, 2, 3, 4, 5])      # like Aspect Ratio
B_values = np.array([1, 2, 3, 4, 5])      # like Mach number

# Create figure
fig, ax = plt.subplots(figsize=(9, 7))

# === Constant A lines (sweep B) - Blue ===
for A in A_values:
    x_list = []
    y_list = []
    for B in B_values:
        # Two simple output functions
        x = A * B                    # x = f(A, B)
        y = A**2 + B                 # y = g(A, B)
        
        x_list.append(x)
        y_list.append(y)
    
    ax.plot(x_list, y_list, 'b-', linewidth=2.2, label=f'A = {A}' if A == A_values[0] else "")
    # Label the line
    mid = len(x_list)//2
    ax.text(x_list[mid], y_list[mid]+0.3, f'A={A}', color='blue', fontsize=11, ha='center')

# === Constant B lines (sweep A) - Red dashed ===
for B in B_values:
    x_list = []
    y_list = []
    for A in A_values:
        x = A * B
        y = A**2 + B
        
        x_list.append(x)
        y_list.append(y)
    
    ax.plot(x_list, y_list, 'r--', linewidth=2, label=f'B = {B}' if B == B_values[0] else "")
    # Label the line
    mid = len(x_list)//2
    ax.text(x_list[mid]+0.3, y_list[mid], f'B={B}', color='red', fontsize=11, fontweight='bold')

# Labels and formatting
ax.set_xlabel('X = A × B', fontsize=12)
ax.set_ylabel('Y = A² + B', fontsize=12)
ax.set_title('Simple Carpet Plot Example\n(Blue = constant A, Red = constant B)', fontsize=14)
ax.grid(True, linestyle='--', alpha=0.6)
ax.legend(title='Curve Families', loc='upper left')

plt.tight_layout()
plt.show()