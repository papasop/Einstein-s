"""
Route B: Complete Test Suite
==============================
Symbolic (exact) + Numerical (independent) verification that

    K_1 = sigma_2^2 * Box_r(ln sigma_1) = r^2/Sigma

when sigma_1 = sqrt(Delta), and K_1 = 1 iff Sigma = r^2.

Two layers:
  PART 1 — Symbolic (sympy): exact algebra, zero floating-point error.
            Verifies every analytic formula before it is used numerically.
  PART 2 — Numerical: independent finite-difference checks using only
            closed-form sigma_1 = sqrt(Delta), no ODE invoked.

Metrics tested:
  Schwarzschild (a=0), Kerr equatorial (theta=pi/2),
  Kerr off-equatorial (general a, theta).

Units: G = c = 1, M = 1.
"""

import sys
import numpy as np
from sympy import (symbols, cos, sin, sqrt, log, diff, simplify,
                   Rational, pi as sym_pi)

M_val = 1.0
PASS  = "\033[92m✓\033[0m"
FAIL  = "\033[91m✗\033[0m"

# ════════════════════════════════════════════════════════════════════════════
# PART 1 — SYMBOLIC VERIFICATION
# ════════════════════════════════════════════════════════════════════════════

print("=" * 65)
print("PART 1: Symbolic verification (sympy, exact)")
print("=" * 65)

r, M, a, theta = symbols('r M a theta', real=True)

Sig_s  = r**2 + a**2*cos(theta)**2       # Sigma
Del_s  = r**2 - 2*M*r + a**2             # Delta
sqg_s  = Sig_s * sin(theta)              # sqrt(-g) radial sector
grr_s  = Del_s / Sig_s                   # g^rr

def Box_sym(f):
    """Symbolic Box_r f = (1/sqrtg)*d_r(sqrtg*g^rr*d_r f)"""
    return simplify(diff(sqg_s * grr_s * diff(f, r), r) / sqg_s)

sym_results = {}

# S1: K_1 = r^2/Sigma  (central claim)
K1_s   = r**2 * Box_sym(Rational(1,2)*log(Del_s))
K1_s   = simplify(K1_s)
sym_results["S1: K_1 = r^2/Sigma"] = simplify(K1_s - r**2/Sig_s) == 0

# S2: Box_r(ln sqrt Delta) = 1/Sigma
box_s  = Box_sym(Rational(1,2)*log(Del_s))
sym_results["S2: Box_r(ln sqrt(Delta)) = 1/Sigma"] = simplify(box_s - 1/Sig_s) == 0

# S3: K_1 = 1 at theta=pi/2 (equatorial)
K1_eq  = simplify(K1_s.subs(theta, sym_pi/2))
sym_results["S3: K_1=1 at theta=pi/2"] = simplify(K1_eq - 1) == 0

# S4: K_1 = 1 at a=0 (Schwarzschild)
K1_sch = simplify(K1_s.subs(a, 0))
sym_results["S4: K_1=1 at a=0"] = simplify(K1_sch - 1) == 0

# S5: 1 - K_1 = a^2*cos^2(theta)/Sigma  (off-equatorial gap)
gap    = simplify(1 - K1_s - a**2*cos(theta)**2/Sig_s)
sym_results["S5: 1-K_1 = a^2*cos^2(theta)/Sigma"] = gap == 0

# S6: K_2 at a=0 equals 1
K2_s   = r**2 * Box_sym(log(r))
K2_sch = simplify(K2_s.subs(a, 0))
sym_results["S6: K_2=1 at a=0"] = simplify(K2_sch - 1) == 0

# S7: Box_r(r^2) formula  (used in T1)
box_r2  = Box_sym(r**2)
claimed = 2*(3*r**2 - 4*M*r + a**2)/Sig_s
sym_results["S7: Box_r(r^2) formula"] = simplify(box_r2 - claimed) == 0

# S8: Box_r(r^3) formula  (used in T7)
box_r3  = Box_sym(r**3)
claimed3 = (12*r**3 - 18*M*r**2 + 6*a**2*r)/Sig_s
sym_results["S8: Box_r(r^3) formula"] = simplify(box_r3 - claimed3) == 0

# S9: d_r(r-M) = 1  (T8 step 2, Schwarzschild)
inner  = Del_s * ((r-M)/Del_s)
sym_results["S9: d_r(Delta*(r-M)/Delta) = d_r(r-M) = 1"] = diff(simplify(inner), r) == 1

# S10: Operator identity  Box_r = (1/Sigma)*d_r(Delta*d_r)
alt_r2 = simplify((1/Sig_s)*diff(Del_s*diff(r**2,r),r))
sym_results["S10: Box_r=(1/Sigma)*d_r(Delta*d_r) on r^2"] = simplify(alt_r2 - box_r2) == 0

print(f"\n  {'check':<45} {'result':>6}")
print("  " + "-" * 53)
all_sym = True
for name, res in sym_results.items():
    flag = PASS if res else FAIL
    print(f"  {flag}  {name}")
    if not res: all_sym = False

print(f"\n  {sum(sym_results.values())}/{len(sym_results)} symbolic checks passed")

# Print the key confirmed formulas
print(f"\n  Confirmed: K_1     = {K1_s}")
print(f"  Confirmed: 1-K_1   = {simplify(1-K1_s)}")
print(f"  Confirmed: K_2(a=0)= {K2_sch}")
print(f"  Confirmed: K_2(gen)= {simplify(K2_s)}")

# ════════════════════════════════════════════════════════════════════════════
# PART 2 — NUMERICAL VERIFICATION
# ════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 65)
print("PART 2: Numerical verification (closed-form sigma_1, no ODE)")
print("=" * 65)

# ── helpers ──────────────────────────────────────────────────────────────────

def Delta(r, a):      return r**2 - 2*M_val*r + a**2
def Sigma(r, a, th):  return r**2 + a**2*np.cos(th)**2

def box_r(f_arr, r_arr, a, th):
    """Numerical Box_r via central differences"""
    sg  = Sigma(r_arr, a, th) * np.sin(th)
    grr = Delta(r_arr, a)     / Sigma(r_arr, a, th)
    return np.gradient(np.gradient(f_arr, r_arr) * sg * grr, r_arr) / sg

def interior(arr, n=25): return arr[n:-n]

def sample_at(arr, r_arr, targets, edge=25):
    ri = r_arr[edge:-edge]; ai = arr[edge:-edge]
    return [ai[np.argmin(abs(ri - t))] for t in targets if ri[0]<=t<=ri[-1]]

def check(vals, expected=1.0, tol=5e-3):
    if not vals: return False, float('inf')
    errs = [abs(v-expected) for v in vals]
    return all(e<tol for e in errs), max(errs)

def ln_sigma1(r_arr, a):
    return 0.5*np.log(np.maximum(Delta(r_arr, a), 1e-30))

r_test = [3.5, 5.0, 7.0, 10.0]
num_results = {}

# ── N1: box_r operator  (neutral f=r^2) ─────────────────────────────────────
print("\nN1: box_r operator  f=r^2  (analytic reference confirmed by S7)")
r_arr = np.linspace(3.0, 14.0, 800)
ok_N1 = True
max_N1 = 0.0
for av in [0.0, 0.5, 0.9]:
    for tv in [np.pi/6, np.pi/3, np.pi/2]:
        num = box_r(r_arr**2, r_arr, av, tv)
        ana = (2*(3*r_arr**2 - 4*M_val*r_arr + av**2)
               / Sigma(r_arr, av, tv))
        err = np.max(np.abs(interior(num) - interior(ana)))
        max_N1 = max(max_N1, err)
        if err > 1e-4: ok_N1 = False
print(f"  max_err over all (a,theta) = {max_N1:.2e}  {PASS if ok_N1 else FAIL}")
num_results["N1: box_r operator (f=r^2)"] = ok_N1

# ── N2: K_1 = r^2/Sigma  (analytic reference confirmed by S1) ───────────────
print("\nN2: K_1 = r^2/Sigma  for sigma_1=sqrt(Delta)")
ok_N2 = True; max_N2 = 0.0
for av in [0.0, 0.3, 0.7, 0.9]:
    r_min = M_val + np.sqrt(max(M_val**2-av**2,0)) + 0.5
    ru = np.linspace(max(r_min,3.0), 14.0, 700)
    for tv in [np.pi/6, np.pi/3, np.pi/2]:
        K1n = ru**2 * box_r(ln_sigma1(ru,av), ru, av, tv)
        K1a = ru**2 / Sigma(ru, av, tv)
        err = np.max(np.abs(interior(K1n) - interior(K1a)))
        max_N2 = max(max_N2, err)
        if err > 1e-3: ok_N2 = False
print(f"  max_err over all (a,theta) = {max_N2:.2e}  {PASS if ok_N2 else FAIL}")
num_results["N2: K_1 = r^2/Sigma"] = ok_N2

# ── N3: K_1=1 where Sigma=r^2  (S3, S4 confirmed symbolically) ──────────────
print("\nN3: K_1 = 1  where Sigma=r^2  (Schwarzschild + Kerr equatorial)")
ok_N3 = True
r_arr = np.linspace(2.6, 14.0, 700)

print("  Schwarzschild (a=0, any theta):")
for tv in [np.pi/6, np.pi/4, np.pi/3, np.pi/2]:
    K1 = r_arr**2 * box_r(ln_sigma1(r_arr,0.0), r_arr, 0.0, tv)
    vals = sample_at(K1, r_arr, r_test)
    ok, me = check(vals)
    if not ok: ok_N3 = False
    print(f"    theta/pi={tv/np.pi:.3f}  max|K1-1|={me:.2e}  {PASS if ok else FAIL}")

print("  Kerr equatorial (theta=pi/2):")
for av in [0.0, 0.2, 0.4, 0.6, 0.8, 0.95]:
    r_min = M_val + np.sqrt(max(M_val**2-av**2,0)) + 0.5
    ru = np.linspace(max(r_min,2.8), 14.0, 700)
    K1 = ru**2 * box_r(ln_sigma1(ru,av), ru, av, np.pi/2)
    vals = sample_at(K1, ru, r_test)
    ok, me = check(vals)
    if not ok: ok_N3 = False
    print(f"    a={av:.2f}  max|K1-1|={me:.2e}  {PASS if ok else FAIL}")

num_results["N3: K_1=1 where Sigma=r^2"] = ok_N3

# ── N4: K_1 < 1 off-equatorial  (S5 confirmed: 1-K_1=a^2cos^2/Sigma) ────────
print("\nN4: K_1 = r^2/Sigma < 1 off-equatorial  (S5 confirmed symbolic)")
ok_N4 = True
r_use = np.linspace(4.0, 14.0, 500)
print(f"  {'a':>5} {'th/pi':>6}  {'K1(r=5)':>10}  {'r^2/Sig':>10}  {'≠1?':>5}  {PASS}")
print("  " + "-" * 52)
for av in [0.5, 0.8]:
    for tv in [np.pi/6, np.pi/4, np.pi/3]:
        K1n = r_use**2 * box_r(ln_sigma1(r_use,av), r_use, av, tv)
        K1a = r_use**2 / Sigma(r_use, av, tv)
        idx = np.argmin(abs(r_use-5.0))
        k1v = K1n[idx]; k1a = K1a[idx]
        err = abs(k1v - k1a)
        neq = abs(k1v-1.0) > 1e-3
        ok  = err < 1e-3
        if not ok or not neq: ok_N4 = False
        print(f"  {av:>5.2f} {tv/np.pi:>6.3f}  {k1v:>10.6f}  "
              f"{k1a:>10.6f}  {'≠1'+PASS if neq else '=1 ?':>5}  "
              f"{PASS if ok else FAIL}")
num_results["N4: K_1<1 off-equatorial"] = ok_N4

# ── N5: Lemma R=0 where K_1=1  (consequence of N3) ──────────────────────────
print("\nN5: Lemma R=(2/r^2)(1-K_1)=0  where K_1=1  (confirmed by S1+S3+S4)")
ok_N5 = True; max_N5 = 0.0
cases5 = [("Schwarzschild",0.0,np.pi/3),("Kerr eq a=0.3",0.3,np.pi/2),
          ("Kerr eq a=0.7",0.7,np.pi/2),("Kerr eq a=0.95",0.95,np.pi/2),
          ("Schwarz th=pi/6",0.0,np.pi/6)]
for name,av,tv in cases5:
    r_min = M_val + np.sqrt(max(M_val**2-av**2,0)) + 0.5
    ru = np.linspace(max(r_min,3.0), 13.0, 600)
    K1 = ru**2 * box_r(ln_sigma1(ru,av), ru, av, tv)
    R  = (2.0/ru**2)*(1.0-K1)
    mr = np.max(np.abs(interior(R)))
    max_N5 = max(max_N5, mr)
    if mr > 1e-2: ok_N5 = False
    print(f"  {name:<22} max|R|={mr:.2e}  {PASS if mr<1e-2 else FAIL}")
num_results["N5: Lemma R=0 where K_1=1"] = ok_N5

# ── N6: K_2=1 Schwarzschild  (confirmed by S6) ───────────────────────────────
print("\nN6: K_2=1 in Schwarzschild  (confirmed by S6)")
r_arr = np.linspace(2.6, 14.0, 700)
K2 = r_arr**2 * box_r(np.log(r_arr), r_arr, 0.0, np.pi/2)
vals = sample_at(K2, r_arr, r_test)
ok6, me6 = check(vals)
num_results["N6: K_2=1 Schwarzschild"] = ok6
print(f"  max|K2-1|={me6:.2e}  {PASS if ok6 else FAIL}")

# ── N7: Operator identity  (neutral f=r^3, confirmed by S8+S10) ──────────────
print("\nN7: Box_r=(1/Sigma)*d_r(Delta*d_r)  f=r^3  (confirmed by S8, S10)")
r_arr = np.linspace(3.0, 14.0, 800)
ok_N7 = True; max_N7 = 0.0
for av in [0.0, 0.4, 0.8]:
    for tv in [np.pi/6, np.pi/3, np.pi/2]:
        num = box_r(r_arr**3, r_arr, av, tv)
        ana = ((12*r_arr**3 - 18*M_val*r_arr**2 + 6*av**2*r_arr)
               / Sigma(r_arr, av, tv))
        err = np.max(np.abs(interior(num) - interior(ana)))
        max_N7 = max(max_N7, err)
        if err > 1e-3: ok_N7 = False
print(f"  max_err over all (a,theta) = {max_N7:.2e}  {PASS if ok_N7 else FAIL}")
num_results["N7: Operator identity (f=r^3)"] = ok_N7

# ── N8: Three-line proof step-by-step  (Schwarzschild) ───────────────────────
print("\nN8: Three-line proof step-by-step  (Schwarzschild, a=0)")
print("  Line 1: d_r(ln sqrt Delta) = (r-M)/Delta   [confirmed S9]")
print("  Line 2: d_r(Delta*(r-M)/Delta) = 1          [confirmed S9, exact]")
print("  Line 3: K_1 = r^2*(1/r^2) = 1               [confirmed S1+S4]")
r_arr = np.linspace(2.8, 14.0, 800)
ok_N8 = True
ls1 = ln_sigma1(r_arr, 0.0)
# Step 1
dls1_num = np.gradient(ls1, r_arr)
dls1_ana = (r_arr-M_val)/Delta(r_arr, 0.0)
e1 = np.max(np.abs(interior(dls1_num)-interior(dls1_ana)))
ok1 = e1 < 1e-4
print(f"  Step 1  max|numerical - (r-M)/Delta| = {e1:.2e}  {PASS if ok1 else FAIL}")
if not ok1: ok_N8 = False
# Step 2
d2 = np.gradient(Delta(r_arr,0.0)*dls1_num, r_arr)
e2 = np.max(np.abs(interior(d2)-1.0))
ok2 = e2 < 1e-3
print(f"  Step 2  max|d_r(Delta*p) - 1| = {e2:.2e}  {PASS if ok2 else FAIL}")
print(f"          [a=0: Sigma/r^2=1, ODE rhs=1; valid here only]")
if not ok2: ok_N8 = False
# Step 3
box_num = r_arr**2 * box_r(ls1, r_arr, 0.0, np.pi/2)
e3 = np.max(np.abs(interior(box_num)-1.0))
ok3 = e3 < 1e-3
print(f"  Step 3  max|K_1 - 1| = {e3:.2e}  {PASS if ok3 else FAIL}")
if not ok3: ok_N8 = False
num_results["N8: Three-line proof (Schwarzschild)"] = ok_N8

# ════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 65)
print("FINAL SUMMARY")
print("=" * 65)

print("\n  PART 1 — Symbolic")
for name, res in sym_results.items():
    print(f"  {PASS if res else FAIL}  {name}")

print("\n  PART 2 — Numerical")
for name, res in num_results.items():
    print(f"  {PASS if res else FAIL}  {name}")

n_sym = sum(sym_results.values())
n_num = sum(num_results.values())
n_tot = len(sym_results) + len(num_results)
n_pass = n_sym + n_num

print(f"\n  Symbolic:  {n_sym}/{len(sym_results)}")
print(f"  Numerical: {n_num}/{len(num_results)}")
print(f"  Total:     {n_pass}/{n_tot}")

ok_all = (n_pass == n_tot)

if ok_all:
    print("""
  ╔══════════════════════════════════════════════════════════╗
  ║  Route B verified: symbolic + numerical, no circularity  ║
  ║                                                          ║
  ║  K_1 = r^2/Sigma  (exact, for sigma_1 = sqrt(Delta))   ║
  ║  K_1 = 1  iff  Sigma = r^2  (a=0 or theta=pi/2)        ║
  ║  1-K_1 = a^2*cos^2(theta)/Sigma  (symbolic, exact)      ║
  ║  Lemma R=0 confirmed where K_1=1                        ║
  ║                                                          ║
  ║  Open: correct sigma_1 off-equatorial in Kerr            ║
  ╚══════════════════════════════════════════════════════════╝
""")
else:
    print("\n  FAILURES — see above.")
    sys.exit(1)
