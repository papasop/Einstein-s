# From Cost to Curvature: Direct Derivation of Einstein's Equations

## Research Roadmap — Y.Y.N. Li, April 2026

---

## 0. The Goal

Derive G_μν = 8πT_μν directly from the information-time cost function,
without passing through Jacobson's thermodynamic bridge.

Current path:   cost → G(2×2) → signature + Clausius → [Jacobson] → Einstein
Target path:    cost[g; δg] → curvature → G_μν = 8πT_μν (one step)

---

## 1. Historical Background: Five Routes to Einstein's Equations

### Route 1: Einstein (1915) — Geometric
- Input: equivalence principle + general covariance + Riemannian geometry
- Method: physical reasoning + trial and error (1911–1915)
- Output: G_μν = 8πT_μν
- Character: geometry tells matter how to move, matter tells geometry how to curve
- Status: the standard derivation

### Route 2: Hilbert (1915) — Variational
- Input: S = ∫(R + L_matter)√(-g) d⁴x
- Method: δS/δg^μν = 0
- Output: G_μν = 8πT_μν
- Character: Einstein's equations as Euler-Lagrange equations
- Advantage: unique — R is the only second-order scalar from the metric
- Limitation: why this action? Why R?

### Route 3: Jacobson (1995) — Thermodynamic
- Input: Lorentzian signature + Clausius δQ = TδS + S = A/4
- Method: apply to all local Rindler horizons
- Output: G_μν = 8πT_μν
- Character: gravity is not fundamental — it's thermodynamic
- Limitation: signature, Clausius, S=A/4 all assumed

### Route 4: Padmanabhan (2010) — Emergent
- Input: horizon thermodynamics + holographic equipartition
- Method: N_sur - N_bulk drives expansion
- Output: Einstein's equations (and beyond, to Lanczos-Lovelock)
- Character: spacetime itself is emergent
- Limitation: still assumes thermodynamic structure

### Route 5: Verlinde (2011) — Entropic
- Input: holographic screens + entropic force F = TΔS/Δx
- Method: derive Newton's law, then extend to GR
- Output: Newton's gravity (GR extension incomplete)
- Character: gravity as entropic force
- Limitation: GR derivation not rigorous; many criticisms

### What all five share
- All require external input about either geometry or thermodynamics
- None derives the metric itself from a more primitive concept
- None explains why Lorentzian signature

### What K=1 adds (current, via Jacobson)
- Derives signature from cost asymmetry (R+E+T)
- Derives Clausius from OU dynamics on {K=1}
- Reduces Jacobson's 3 external inputs to 1 (S ∝ A)
- But still passes through the thermodynamic bridge

### What the direct route would add
- Derive curvature directly from cost structure
- Define matter from cost deviations
- Einstein's equations as cost-functional stationarity
- No thermodynamic intermediary

---

## 2. The Conceptual Framework

### 2.1 From point-cost to field-cost

Current K=1:
  d(x; δx) — cost of displacing from x by δx
  G_ij = ∂²(dt_info²)/∂δx^i ∂δx^j — Hessian at a point
  → 2×2 matrix at each Rindler horizon

Extension needed:
  D[g; δg] — cost of deforming metric g by δg
  G_{μναβ} = δ²D/δ(δg^μν)δ(δg^αβ) — Hessian on metric space
  → DeWitt supermetric on the space of 4-metrics

Key insight:
  The DeWitt supermetric G^{μναβ} already exists in canonical GR:
    G^{μναβ} = (1/2)(g^μα g^νβ + g^μβ g^να) - g^μν g^αβ
  It has Lorentzian signature on superspace (one negative direction = conformal mode)
  
  Question: can this be derived from a cost function on metric space?

### 2.2 Three-level structure

Level 0 (current): cost on R² → metric G (2×2)
Level 1 (needed):  cost on spacetime → metric g_μν (4×4)
Level 2 (needed):  cost on metric space → curvature R_μν

The direct derivation requires all three levels.

### 2.3 The key hypothesis

If K(x) = x^T G x measures temporal cost at a point,
then the field-theory analogue is:

  K[g] = ∫ g^μν R_μν √(-g) d⁴x = ∫ R √(-g) d⁴x

This is the Einstein-Hilbert action!

Hypothesis: R is the field-theory version of K = x^T G x.

  K = x^T G x     (point):  metric acting on coordinates
  R = g^μν R_μν   (field):  metric acting on its own curvature

Both are "the metric measuring itself":
  K: metric measures the cost of a displacement
  R: metric measures its own curvature

If K=1 at the point level → self-consistent metric
then R = const at the field level → self-consistent spacetime?

This would give: δ∫R√(-g)d⁴x = 0 → Einstein's equations (vacuum)

---

## 3. The Roadmap: Four Steps

### Step 1: Cost on a manifold (promote G to g_μν)

Task: let cost depend on position
  Current: G = const (one Rindler horizon)
  Extension: G → g_μν(x) (varies across spacetime)

Concrete: 
  Take G(r) = diag(-σ₁(r)², 1) where σ₁ = σ₁(r)
  The cost d² = δx^T G(r) δx now depends on where you are
  K(x; r) = δx^T G(r) δx

What this gives:
  A family of 2×2 metrics parametrized by position
  → effectively a metric on a 2D spacetime
  → Christoffel symbols from ∂G/∂r
  → Riemann tensor from ∂²G/∂r²

Verification test:
  For Schwarzschild: σ₁(r) = r√(1-2M/r)
  Does the cost-derived curvature match R = 0 (vacuum)?
  → One calculation. Can be done immediately.

### Step 2: Curvature from cost derivatives

Task: express Riemann tensor as cost derivatives

The hierarchy:
  0th order: d(x; δx) → cost function
  2nd order: ∂²d²/∂δx² → metric g_μν (Hessian)
  3rd order: ∂g/∂x → Christoffel Γ^ρ_μν
  4th order: ∂Γ/∂x → Riemann R^ρ_σμν

Question: can R^ρ_σμν be expressed directly in terms of cost?

If d²(x; δx) = g_μν(x) δx^μ δx^ν (to leading order), then:
  g_μν = ∂²d²/∂δx^μ ∂δx^ν |_{δx=0}
  Γ^ρ_μν = (1/2)g^ρσ(∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
  R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ

This is just standard differential geometry applied to the cost-derived metric.
Nothing new mathematically — but the interpretation is new:
  "Curvature = how the cost structure varies across spacetime"

### Step 3: K=1 field equation

Task: promote K=1 from a point condition to a field equation

Point level:
  K = x^T G x = 1
  V = (1/2)(K-1)² → restoring dynamics

Field level candidates:
  
  (a) R = 0 (vacuum Einstein)
      "The scalar curvature vanishes" = "spacetime is self-consistent"
      Analogy: K=1 means "cost = 1"; R=0 means "total curvature cost = 0"
      Problem: only vacuum. No matter.

  (b) R = const (cosmological constant)
      "Curvature is uniform" = "self-consistency allows a constant offset"
      K=1 → R = Λ? Λ from the "1" in K=1?

  (c) δ∫(R-2Λ)√(-g)d⁴x = 0
      Variation → G_μν + Λg_μν = 0 (vacuum with Λ)
      The "self-consistent" spacetimes are those that extremize
      the total curvature cost.

  (d) With matter: K[g] ≠ 1 where matter is present
      K-1 ≠ 0 → "cost deviation" = matter
      → T_μν ∝ (K[g]-1) in some sense
      → G_μν = 8π × (deviation from self-consistency)

Option (d) is the most promising:
  T_μν = "the metric's failure to be self-consistent"
  Einstein's equations = "curvature responds to self-inconsistency"

### Step 4: Matter from cost deviation

Task: define T_μν from the cost function

In the point-level K=1:
  K = 1: vacuum (self-consistent)
  K ≠ 1: deviation → generates restoring force → V = (1/2)(K-1)²

Field-level analogue:
  If "vacuum" = "metric is self-consistent" (R_μν = 0 for Ricci-flat)
  then "matter" = "metric is NOT self-consistent" (R_μν ≠ 0)
  
  T_μν := -(1/8π)(R_μν - (1/2)Rg_μν)
  
  This is just the definition of T_μν from Einstein's equations!
  But the interpretation changes:
    Standard GR: T_μν is independently defined (from matter Lagrangian)
    Cost framework: T_μν IS the self-inconsistency of the metric

  This is radical but not new — Sakharov (1967) proposed "induced gravity"
  where the metric is fundamental and matter is emergent.

---

## 4. Verification Plan

### Test 1: Schwarzschild (immediate)
  σ₁(r) = r√(1-2M/r)
  G(r) = diag(-σ₁(r)², 1) in (τ,ℓ) coordinates
  Compute Christoffel → Riemann → Ricci
  Check: R_μν = 0 (vacuum)?
  → This should work — it's just standard GR in different notation
  → But verifies that cost-derived metric gives correct curvature

### Test 2: Symplectic curvature formula
  From your symplectic gravity paper:
    R = -2□lnσ₁ + 2/r²
  Check: does this match the cost-derived Riemann tensor?
  → Connects K=1 to your existing results

### Test 3: FRW cosmology
  σ₁(t) = a(t) (scale factor)
  G(t) = diag(-a(t)², 1)
  K[a] = ?
  Does K=1 at field level → Friedmann equation?
  → This would be a strong test

### Test 4: Perturbation theory
  g_μν = η_μν + h_μν (weak field)
  Cost of h_μν perturbation?
  Linearized Einstein equations from cost variation?
  → Most tractable analytically

---

## 5. Key References

- Einstein (1915): field equations from physical reasoning
- Hilbert (1915): field equations from variational principle
- Sakharov (1967): induced gravity — geometry is fundamental, matter is emergent
- Jacobson (1995): Einstein from thermodynamics
- Padmanabhan (2010): emergent spacetime program
- Verlinde (2011): entropic gravity
- DeWitt (1967): supermetric on the space of 3-geometries
- Amari (2016): information geometry (Riemannian Fisher metric)
- Li (2026): K=1 Chronogeometrodynamics (the present framework)

---

## 6. Risk Assessment

### What could go wrong

1. Cost-derived curvature may be trivially equivalent to standard GR
   → Just renaming, no new content
   → Then the direct route adds nothing over Jacobson route

2. K=1 at field level may not give Einstein's equations
   → May give a different field equation
   → Could be interesting (modified gravity) or wrong

3. Matter definition from cost deviation may be circular
   → "T_μν = G_μν/8π" is just the definition
   → Need independent content

4. 4D extension of K=1 may require arbitrary choices
   → Symplectic pairing ambiguity (dim Skew₄ = 6)
   → May lose uniqueness

### What could go right

1. R = field-level K is not just analogy but derivable
   → Cost functional's Hessian → DeWitt supermetric
   → Second variation → Riemann tensor
   → This would be genuinely new

2. Matter emerges naturally as cost deviation
   → T_μν has information-theoretic meaning
   → "Matter = metric's self-inconsistency"
   → Connects to Wheeler's "matter from geometry"

3. Cosmological constant appears as the "1" in K=1
   → Λ = the field-level normalisation constant
   → Would give Λ a geometric origin
   → Major open problem in physics

---

## 7. Immediate Next Steps

### Step 1 (today): Schwarzschild curvature test
  Compute curvature from position-dependent cost
  Verify R=0 for vacuum Schwarzschild

### Step 2 (this week): linearized theory
  Perturbation h_μν around flat space
  Cost of perturbation → linearized Einstein?

### Step 3 (if Step 2 succeeds): field-level K
  Define K[g] = ∫R√(-g)d⁴x / ∫√(-g)d⁴x (averaged curvature)
  Study K=1 as field equation
  Compare with Einstein-Hilbert variation

### Step 4 (if Step 3 succeeds): write paper
  "Einstein's Equations from Information-Time Cost"
  → Direct derivation, no Jacobson
  → Submit to PRL (if successful, this is PRL-level)

---

## 8. The Central Question

At the point level:
  K = x^T G x — the metric measuring the cost of a displacement

At the field level:
  R = g^μν R_μν — the metric measuring its own curvature

Are these the same operation at different scales?

If yes: Einstein's equations are the field-theoretic version of K=1.
If no: the analogy is formal, and Jacobson remains the only bridge.

One calculation will tell.
