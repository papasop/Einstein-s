"""
Derivation of K_angular from cost principle
============================================================
Key insight: K_angular = σ₂² □ln σ₂ 
(same structure as K_field = σ₂² □ln σ₁)
"""
import sympy as sp
from sympy import *

r, M = symbols('r M', positive=True)
f = 1 - 2*M/r

print("=" * 70)
print("CLAIM: K_angular = σ₂² □ln σ₂")
print("=" * 70)

# σ₂ = r
sigma2 = r
ln_sigma2 = ln(sigma2)  # = ln(r)

# □φ = (1/r²) d/dr [r²f · dφ/dr]  for static spherically symmetric
d_ln_s2 = diff(ln_sigma2, r)  # = 1/r
print(f"\n  σ₂ = r")
print(f"  ln σ₂ = ln r")
print(f"  d(ln σ₂)/dr = {d_ln_s2}")

inner = r**2 * f * d_ln_s2  # = r²f/r = rf
inner_simplified = simplify(inner)
print(f"  r²f·(ln σ₂)' = {inner_simplified}")

d_inner = diff(inner, r)
d_inner_simplified = simplify(d_inner)
print(f"  d/dr[r²f·(ln σ₂)'] = {d_inner_simplified}")

box_ln_s2 = simplify(d_inner / r**2)
print(f"  □ln σ₂ = {box_ln_s2}")

K_angular_from_box = simplify(r**2 * box_ln_s2)
print(f"\n  σ₂²·□ln σ₂ = r²·□ln σ₂ = {K_angular_from_box}")

# Compare with rf' + f
K_angular_standard = simplify(r * diff(f, r) + f)
print(f"  rf' + f = {K_angular_standard}")

match = simplify(K_angular_from_box - K_angular_standard) == 0
print(f"\n  σ₂²·□ln σ₂ = rf' + f ?  {match}")

if match:
    print("  ✓ K_angular = σ₂²·□ln σ₂  CONFIRMED")
else:
    print(f"  ✗ Difference: {simplify(K_angular_from_box - K_angular_standard)}")

print("\n" + "=" * 70)
print("VERIFICATION: K_angular = 1 for Schwarzschild")
print("=" * 70)

K_ang_schw = simplify(K_angular_from_box)
print(f"\n  K_angular = σ₂²·□ln σ₂ = {K_ang_schw}")
print(f"  K_angular = 1 ?  {K_ang_schw == 1}")

print("\n" + "=" * 70)
print("THE UNIFIED PRINCIPLE")
print("=" * 70)

# K_field = σ₂² □ln σ₁
sigma1 = r * sqrt(f)
ln_sigma1 = ln(r) + Rational(1,2)*ln(f)
d_ln_s1 = diff(ln_sigma1, r)
inner1 = r**2 * f * d_ln_s1
K_field = simplify(diff(inner1, r))

print(f"""
  σ₁ = r√f,  σ₂ = r

  K_field   = σ₂²·□ln σ₁ = {K_field}
  K_angular = σ₂²·□ln σ₂ = {K_angular_from_box}

  UNIFIED PRINCIPLE:
    σ₂²·□ln σᵢ = 1  for each symplectic eigenvalue σᵢ

    i=1: σ₂²·□ln σ₁ = 1  →  R = 0
    i=2: σ₂²·□ln σ₂ = 1  →  R_θθ = 0

  Together: R_μν = 0 (vacuum Einstein equations)
""")

# Verify for general f
print("=" * 70)
print("GENERAL f(r): verify K_angular = σ₂²·□ln σ₂")
print("=" * 70)

f_gen = Function('f')
ln_s2_gen = ln(r)
d_ln_s2_gen = diff(ln_s2_gen, r)
inner_gen = r**2 * f_gen(r) * d_ln_s2_gen
K_ang_gen = simplify(diff(inner_gen, r))
K_ang_standard_gen = r*diff(f_gen(r), r) + f_gen(r)

print(f"\n  σ₂²·□ln σ₂ = {K_ang_gen}")
print(f"  rf' + f     = {K_ang_standard_gen}")
print(f"  Equal?  {simplify(K_ang_gen - K_ang_standard_gen) == 0}")

print("\n" + "=" * 70)
print("WHAT THIS MEANS")
print("=" * 70)

print("""
  BEFORE (old paper):
    K_field   = σ₂²·□ln σ₁ = 1  ← derived from cost
    K_angular = rf' + f = 1      ← DEFINED to match R_θθ=0 (NOT derived)
    
    Gap: K_angular was imposed, not derived.

  AFTER (this derivation):
    K_field   = σ₂²·□ln σ₁ = 1  ← instance of unified principle
    K_angular = σ₂²·□ln σ₂ = 1  ← SAME principle, different eigenvalue
    
    No gap: both are instances of
      "σ₂²·□ln σ = 1 for each symplectic eigenvalue σ"

  The principle:
    □ln σ = 1/σ₂²  for every eigenvalue σ

  Physical meaning:
    "Every symplectic eigenvalue propagates at the angular curvature rate."
    "The d'Alembertian of each log-eigenvalue equals the angular cost 1/r²."
    
  This is cost self-consistency applied to EACH eigenvalue separately,
  not just to σ₁.
  
  K_angular was ALWAYS σ₂²□lnσ₂ — we just didn't recognize it.
""")

# Final: verify the equivalence theorem with the new derivation
print("=" * 70)
print("EQUIVALENCE THEOREM (now fully derived)")
print("=" * 70)

C1, C2 = symbols('C1 C2')
f_general = 1 - C1/r - C2/r**2

# K_field for general solution
ln_s1_gen2 = ln(r) + Rational(1,2)*ln(f_general)
d1_gen2 = diff(ln_s1_gen2, r)
inner1_gen2 = r**2 * f_general * d1_gen2
K_f_gen2 = simplify(diff(inner1_gen2, r))

# K_angular for general solution  
inner2_gen2 = r**2 * f_general / r  # = r·f
K_a_gen2 = simplify(diff(inner2_gen2, r) * r**2 / r**2)
# Actually just compute directly
K_a_gen2 = simplify(r*diff(f_general, r) + f_general)

print(f"\n  f = 1 - C₁/r - C₂/r²")
print(f"  K_field   = σ₂²·□ln σ₁ = {K_f_gen2}")
print(f"  K_angular = σ₂²·□ln σ₂ = {K_a_gen2}")
print(f"\n  K_field = 1: always ✓ (for this family)")
print(f"  K_angular = 1: requires C₂ = {sp.solve(K_a_gen2 - 1, C2)}")
print(f"\n  Unified principle σ₂²□lnσ = 1 for both eigenvalues")
print(f"  → f = 1 - 2M/r (Schwarzschild, unique)")
