import math

# main wing
Sw = 543.03796
bw = 35.00000
MACw = 17.06951   

# horizontal tail 
Sh = 77.27030
bh = 19.36488
MACh = 4.18160

# vertical tails 
Sv_total = 86.21542
bv_span = 16.00000       
bv_proj = 15.76653      
MACv = 5.78128

# configuration
n_stabs = 2  
n_fins  = 2   

lh = 15.41  # ft
lv = 12.88  # ft

# control surfaces
rudder_chord_frac = 0.30
rudder_span_frac  = 0.80
elev_chord_frac   = 0.30
elev_span_frac    = 0.80

ARw_calc = bw**2 / Sw
ARh_calc = bh**2 / Sh

# vertical tail AR 
ARv_calc = bv_proj**2 / Sv_total

cbar_w = Sw / bw
cbar_h = Sh / bh

# vertical tail mean chord
cbar_v_span = Sv_total / bv_span
cbar_v_proj = Sv_total / bv_proj

# surface areas
Sh_each = Sh / n_stabs
Sv_each = Sv_total / n_fins

# elevators
elev_chord = elev_chord_frac * cbar_h
elev_span  = elev_span_frac  * (bh / n_stabs)  

# rudders
rudder_chord = rudder_chord_frac * cbar_v_span
rudder_height = rudder_span_frac * bv_proj

# volume coefficients
def volume_coeffs(Sw, bw, MACw, Sh, lh, Sv_total, lv):
    Vh = (lh * Sh) / (MACw * Sw)
    Vv = (lv * Sv_total) / (bw * Sw)
    return Vh, Vv

print("Main Wing:")
print(f"  Sw = {Sw:.5f} ft^2, bw = {bw:.5f} ft, MACw = {MACw:.5f} ft")
print(f"  AR = {ARw_calc:.5f}")
print(f"  mean chord Sw/bw = {cbar_w:.5f} ft\n")

print("Horizontal Tail:")
print(f"  Sh = {Sh:.5f} ft^2, bh = {bh:.5f} ft, MACh = {MACh:.5f} ft")
print(f"  AR = {ARh_calc:.5f}")
print(f"  mean chord Sh/bh = {cbar_h:.5f} ft")
print(f"  per stabilizer area = {Sh_each:.5f} ft^2\n")

print("Vertical Tails:")
print(f"  Sv_total = {Sv_total:.5f} ft^2, bv_span = {bv_span:.5f} ft, bv_proj = {bv_proj:.5f} ft")
print(f"  AR = {ARv_calc:.5f}")
print(f"  mean chord (Sv/span) = {cbar_v_span:.5f} ft")
print(f"  mean chord (Sv/proj span) = {cbar_v_proj:.5f} ft")
print(f"  per fin area = {Sv_each:.5f} ft^2\n")

print(f"Elevator: ~{elev_chord_frac*100:.0f}% chord -> {elev_chord:.3f} ft, "
      f"~{elev_span_frac*100:.0f}% span (per stab) -> {elev_span:.3f} ft")
print(f"Rudder:  ~{rudder_chord_frac*100:.0f}% chord -> {rudder_chord:.3f} ft, "
      f"~{rudder_span_frac*100:.0f}% height -> {rudder_height:.3f} ft\n")

Vh, Vv = volume_coeffs(Sw, bw, MACw, Sh, lh, Sv_total, lv)
print("Tail volume coefficients:")
print(f"Vh = {Vh:.5f}")
print(f"Vv = {Vv:.5f}")