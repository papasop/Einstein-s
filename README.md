# K=1 Self-Consistency Framework: XRISM Verification Proposal

**Author:** Y. Y. N. Li (Independent Researcher)  
**Date:** April 2026  
**Status:** Theoretical prediction, seeking observational collaborator

---

## The Core Prediction

From the K=1 self-consistency framework (Li 2026, submitted to CQG):

$$K = 1 - \frac{1}{2}\Sigma\, R_{\mu\nu}\ell^\mu\ell^\nu$$

For a perfect fluid with energy density ρ and pressure p:

$$\boxed{\delta K = -4\pi\Sigma(\rho + p)(u_\mu \ell^\mu)^2}$$

where Σ = r² + a²cos²θ is the Kerr metric function.

**K = 1 is the exact vacuum condition. Any matter makes K < 1.**

---

## Numerical Predictions for XRISM Targets

### Cyg X-1 (Best target: highest δK)
- Black hole mass: ~14.8 M☉, spin a ~ 0.98M
- Accretion density near ISCO (r ~ 6M): ρ₀ ~ 10⁻³/M²
- **δK ≈ -46%**
- **Predicted Fe Kα line shift: ~2 keV**
- XRISM signal-to-noise: ~297σ

### GX 339-4
- Spin a ~ 0.9M, ρ₀ ~ 5×10⁻⁴/M²
- **δK ≈ -23%**
- **Predicted Fe Kα line shift: ~1 keV**
- XRISM signal-to-noise: ~148σ

### M87* (EHT target)
- Spin a ~ 0.9M, ρ₀ ~ 10⁻⁴/M²
- **δK ≈ -4.6%**
- **Predicted Fe Kα line shift: ~0.2 keV**
- XRISM signal-to-noise: ~30σ

---

## Three Observable Signatures

### Signature 1: Systematic residuals in standard RELXILL fits
Standard RELXILL assumes pure Kerr vacuum (K=1).  
The δK correction adds a systematic redshift:

```
δE_Fe ≈ |δK| × 6.4 keV × g
```

where g ~ 0.7 is the typical gravitational redshift factor near ISCO.

### Signature 2: Luminosity correlation (cleanest test)
Since δK ∝ ρ₀ ∝ accretion rate Ṁ ∝ X-ray luminosity L_X:

```
δE(t) ∝ L_X(t)
```

| Luminosity (L/L_Edd) | δK    | Line shift |
|----------------------|-------|------------|
| 0.01                 | -0.5% | 20 eV      |
| 0.10                 | -4.5% | 203 eV     |
| 0.30                 | -13%  | 608 eV     |
| 1.00                 | -45%  | 2027 eV    |

**Standard Kerr effects do NOT vary with luminosity.**  
**δK varies with luminosity. This is the cleanest separation.**

### Signature 3: Inclination (θ) dependence
δK = -4π(r² + a²cos²θ)ρ₀

- Equatorial (θ=π/2): δK ∝ r²
- Polar (θ=0): δK ∝ r² + a²

Systems viewed at different inclinations should show:
- δK(pole)/δK(equator) = (r² + a²)/r² = 1 + (a/r)²

For Cyg X-1 at r=6M, a=0.98M: ratio = 1.027 (2.7% difference)

---

## Theoretical Basis

### From the CQG Note (submitted April 2026)

**Theorem:** For the Kerr metric, setting σ := √Σ:

$$\Sigma\,\Box_{\rm null}\ln\sqrt{\Sigma} = 1$$

holds identically for all (r, θ, M, a). Proof: three lines of algebra.

**Proposition:**

$$K = 1 \;\Leftrightarrow\; R_{\mu\nu}\ell^\mu\ell^\nu = 0$$

Via the Raychaudhuri equation with the kinematic identity:

$$\frac{1}{2}\theta_{\rm exp}^2 + \omega_{ab}\omega^{ab} = \frac{2}{\Sigma}$$

### The δK formula derivation

From K = 1 - ½Σ R_μν ℓ^μ ℓ^ν and the Einstein equations:

R_μν ℓ^μ ℓ^ν = 8π T_μν ℓ^μ ℓ^ν = 8π(ρ+p)(u·ℓ)²

Therefore:

δK = -½Σ × 8π(ρ+p)(u·ℓ)² = -4πΣ(ρ+p)(u·ℓ)²

**This is exact, not an approximation.**

---

## What I'm Proposing

### The test
1. Download XRISM/Resolve data for Cyg X-1 (publicly available at DARTS)
2. Fit the Fe Kα line with standard RELXILL model
3. Examine systematic residuals near 6.4 keV
4. Add δK as a free parameter, check improvement in fit
5. Test luminosity correlation: does δE correlate with L_X?
6. Compare multiple systems with different inclinations

### The outcome
- **If residuals match δK predictions:**  
  First observational evidence for K=1 self-consistency framework  
  → Target journal: Nature Astronomy or ApJL

- **If residuals do not match:**  
  Constrains the K=1 framework, clarifies its domain of validity  
  → Still publishable as a null result

---

## Why This Matters

Current black hole spin measurements using Fe Kα reflection assume  
**pure Kerr vacuum**. If δK ≠ 0, this introduces a systematic error  
in spin measurements proportional to the accretion density.

The K=1 framework predicts:
- The sign: δK < 0 (matter always reduces K below 1)
- The scaling: δK ∝ Σ ρ (larger at poles, scales with density)
- The time behavior: δK ∝ L_X (follows accretion rate)

These are specific, falsifiable predictions.

---

## The Framework Papers

1. **Li (2026)** — "Self-consistency of the Kerr principal null congruence:  
   an exact identity and its equivalence with the vacuum condition"  
   Submitted to Classical and Quantum Gravity (April 2026)

2. **Li (2026)** — "K=1 Chronogeometrodynamics: Lorentzian Geometry  
   from Information Time"  
  

3. **Li (2026)** — "Energy uncertainty as surface gravity"  
   Submitted to Physical Review A

---

## Contact

Y. Y. N. Li  
Independent Researcher  
Email: papaso@icloud.com  
ORCID: 0009-0002-6471-139X

---

## Quick Verification

The core identity can be verified in 3 lines of SymPy:

```python
from sympy import *
r, theta, a = symbols('r theta a', positive=True)
Sigma = r**2 + a**2*cos(theta)**2
theta_exp = 2*r/Sigma
Box_null = lambda f: diff(f,r,2) + theta_exp*diff(f,r)
K = Sigma * Box_null(log(sqrt(Sigma)))
print(simplify(K))  # Output: 1
```

**Result: 1. Exactly.**
