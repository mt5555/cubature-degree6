# Degree 6: algebraic cross-check by the critical-point method

This is an **independent cross-check** of the degree-6 theorem in the main README: no 10-point
degree-6 rule has all nodes in the closed triangle. It uses exact algebraic solving instead of an
SOS certificate. **Status: computed with msolve and Singular, but not formally certified.** The
rigorous proof is the SOS certificate in the main README; see "What is not certified" below.

## Setup

* Unknowns: the 8 degree-7 moments y. Equations: the 6 flatness quadrics h(y) = 0, from
  `../cf2.py`. Their solution set V is the set of all 10-point degree-6 rules.
  (`gen6.py`)
* Inside condition: K = {G_g(y) ⪰ 0 for g in {x, y, 1−x−y}}, the three 4×4 localizing LMIs. On V,
  det G_g = (positive factor) · Π wᵢ g(xᵢ). So the boundary of K is where some node lies on an edge line.
* **Structure of V** (Singular, `V.sing`): Krull dimension 2 (a surface), degree 26. The Bézout bound
  is 64; the difference is accounted for by solutions at infinity.

## Argument

V ∩ K is compact. If it were nonempty, a generic linear function ℓ (here
ℓ = 3y₀ − 5y₁ + 7y₂ + 2y₃ − 4y₄ + 6y₅ − y₆ + 8y₇) would attain its maximum on it at some point p.
Depending on which boundaries of K are active at p, p is one of the following:

1. a singular point of V, or a critical point of ℓ on V (no boundary active);
2. a singular point of a boundary curve V ∩ {det G_g = 0}, or a critical point of ℓ on that curve
   (one boundary active);
3. a point of V ∩ {det G_g = 0} ∩ {det G_g′ = 0} (two boundaries active).

All of these sets are finite. Each was computed, and every real point was tested against K:

| Candidate type | System | Complex | Real | Real points in K |
|---|---|---|---|---|
| Singular points of V | h + 6×6 minors of ∂h (Singular, `singV.sing`) | 0 (empty) | 0 | none |
| Critical points of ℓ on V | Lagrange: ∇ℓ = Σλ∇h, h = 0; 14 eqs / 14 unknowns (`lagrange.py`) | 114 | 8 | none; worst whitened eigenvalue ≤ −1.38 |
| Critical points on boundary curves | ∇ℓ = Σλ∇h + μ∇det_g, h = 0, det_g = 0; 15 / 15 (`lagrange_bd.py`) | 90 per edge | 16, 14, 18 | none; each has a clearly negative eigenvalue on its own edge, or fails another edge by ≥ 0.44 |
| Singular points of boundary curves | ∇det_g = Σc∇h, h = 0, det_g = 0 (`sing_bd.py`) | 12 per edge | 4 per edge | none; ≤ −2.4, or far outside the bounding box |
| Corners (two edges active) | h, det_g, det_g′ (Singular, `critB.sing`) | 38 per pair | 12 per pair | none; ≤ −1.8 |

No candidate lies in K, so V ∩ K = ∅. This agrees with the SOS proof, and with large margins.
(Whitened eigenvalues use the scaling where Lebesgue measure gives 1, and K requires ≥ 0.)
`rur_points.py` evaluates the real points from msolve's rational parametrizations at 300 digits
and runs the K test.

## What is not certified

The cross-check takes msolve's output as correct and complete. A formal version would need,
for each of the 10 systems, the same three steps used in `../degree8/`:

1. exact verification of each rational parametrization by substitution;
2. an exact Gröbner-basis count over ℚ (plain Singular `std`) to show no solutions are missing;
3. exact interval checks that each real candidate violates an LMI.

These steps were deliberately skipped, because the SOS proof is already rigorous and was
independently reviewed. The two non-msolve Singular results (V's dimension and degree, and the
corner counts) use `modStd`, which is not certified for upper bounds either.

## Files

* Data: `gen6.py` → `data6.pkl`.
* Systems: `lagrange.py`, `lagrange_bd.py`, `sing_bd.py`, and `critB.sing` (corners), with inputs
  `*.ms` and msolve outputs `*.out`; `V.sing`, `singV.sing`.
* Evaluation: `rur_points.py` (uses `parse_rur.py`).
