# Degree 8: every 15-point rule on the triangle has nodes outside

**Theorem.** On T = {x ≥ 0, y ≥ 0, x + y ≤ 1}, there are exactly **two** 15-point cubature
rules that integrate every polynomial of degree ≤ 8 exactly. The two are mirror images.
Each has 12 nodes inside T and 3 nodes far outside, each outside node carrying weight
1.6×10⁻¹⁰. So **no 15-point degree-8 rule has all nodes in the closed triangle.** The best
known rule with all nodes inside has 16 points (Wandzura–Xiao). Since 15 = dim P₄ is the
lower bound, the minimal number of nodes for a degree-8 rule with all nodes inside is 16.

This is a computer-assisted proof. Unlike degree 6, it needs no SOS certificate: the solution
set is finite, so it is enumerated exactly.

## 1. Reduction to 10 equations in 10 unknowns

This is the degree-6 argument one size up; see the main README, Section 1.

* **Weights.** Any 15-point rule exact to degree 8 satisfies Vᵀ diag(w) V = M₄, the Lebesgue
  moment matrix on P₄ (15×15, positive definite). So the nodes are distinct and P₄-unisolvent,
  and every weight is positive, for any real weights.
* **Flatness.** M₅ has rank 15 = rank M₄. The unknowns are the 10 degree-9 moments
  yₖ = μ(x^(9−k) yᵏ). B = M[P₄, degree-5 monomials] is 15×6, and C = Bᵀ M₄⁻¹ B must be
  Hankel. That gives **10 quadratic equations h(y) = 0.** Node space agrees: 45 unknowns and
  dim P₈ = 45 equations.
* **Nodes inside.** The localizing matrices M₄(g·μ) must be PSD for g in {x, y, 1−x−y}.
  Taking the Schur complement on the fixed positive-definite P₃ block gives three 5×5
  linear matrix inequalities G_g(y) ⪰ 0.
* **Converse.** By Curto–Fialkow, every real solution of h = 0 is a genuine 15-point rule
  with positive weights.

## 2. Exact enumeration (msolve 0.6.5)

* msolve computes a rational parametrization of all complex solutions of h = 0: an
  eliminating polynomial f(t) of degree 16, and yᵢ = −nᵢ(t)/(cᵢ q(t)).
* `verify_rur_exact.py` checks exactly, in rational arithmetic:
  * f is squarefree and gcd(f, q) = 1, so the parametrization gives 16 distinct points;
  * every one of the 10 equations vanishes on all 16 points, i.e. Eₖ(y(t))·q(t)² ≡ 0 mod f;
  * f has exactly 2 real roots (Sturm), so exactly 2 of the 16 solutions are real.
* **Completeness: there are no other solutions.** Singular's plain `std` (`std_Q.sing`) computes
  a Gröbner basis of the ideal over ℚ in exact rational arithmetic. It uses Buchberger's
  algorithm with no modular or probabilistic steps, and takes about 15 minutes. It reports
  Krull dimension 0 and dim_ℚ ℚ[y]/I = 16 (`std_Q.log`). So the ideal has at most 16 solutions
  counted with multiplicity. The 16 verified distinct points are therefore all of them, and
  each has multiplicity 1. Singular's multi-modular `modStd` gives the same answer in 6 s
  (`modstd_Q.log`).
* `real_outside.py`: for each real root, it takes a rational isolating interval for t,
  encloses y(t) in exact interval arithmetic (width ≤ 2×10⁻¹⁰), and finds a rational vector v
  with vᵀ G_g(y) v ≈ −7.6×10⁻⁶ < 0 over the whole enclosure. So each real solution violates a
  localizing condition, meaning some node lies outside T.
* `verify_data8.py` re-derives the problem data independently from the direct 15×15 moment
  matrices:
  * M₄ is positive definite;
  * the msolve input equations are exactly nonzero rational multiples of the Hankel conditions;
  * the Schur complements of the 15×15 localizing matrices equal G_g, and each P₃ block A_g is positive definite.

## 3. The two real rules

| | Rule A | Rule B (mirror image of A) |
|---|---|---|
| nodes inside T | 12 (min edge distance 0.0246, min weight 0.0265) | same |
| nodes outside T | 3: (3.623, −2.969), (−2.969, 0.346), (0.346, 3.623) | mirrored |
| weight of each outside node | 1.606×10⁻¹⁰ | same |
| symmetry | 3-fold rotation, no reflection | same |

Both refine to node-space residual ~10⁻¹¹ with weights summing to 1/2. They are presumably
the "15†" entry in Cools' encyclopedia.

## 4. Numerical cross-checks and a tooling note

* A Newton search in the moment unknowns (4,500 starts) found exactly these two real solutions.
* Node-space least squares with all 15 nodes constrained inside T (2×150 starts): the best
  residual is 0.0947 out of 0.707, far from an exact rule.
* **msolve pitfall.** In prime-field mode, msolve garbles input coefficients larger than the
  prime unless they are reduced mod p first. This gave a spurious degree of 974. With
  coefficients reduced mod p, three primes all give 16, agreeing with the exact computation over ℚ.

## 5. Literature

* Papanicolopulos (2015, 2016) and Witherden–Vincent (2015) prove minimality among *fully
  symmetric* rules (up to degree 14), and report new symmetric rules. None reports a 15-point
  degree-8 rule with all nodes inside.
* Yuan Xu, *Minimal Cubature Rules: Theory and Practice* (Cambridge, 2025), has a chapter on
  the triangle. It has not been checked here; its abstract emphasizes invariant (symmetric) rules.
* The result here covers *all* 15-point rules, with no symmetry assumed.

## Files

`cf8.py` (exact formulation), `verify_data8.py`, `write_msolve.py` (writes `flat8_Q.ms`),
`rur_Q.out` (msolve parametrization), `parse_rur.py`, `verify_rur_exact.py`, `real_outside.py`,
`ysearch8.py`, `analyze8.py`, `nodesearch8.py`, `common8.py`, `reduce_mod.py`, `std_Q.sing`/`std_Q.log`,
`modstd_Q.sing`/`modstd_Q.log`, and logs `verify_data8.log`, `verify_rur_exact.log`, `real_outside.log`.

**Reproducing** (needs `apt install msolve singular`, plus python with sympy and numpy):
`python3 verify_data8.py`, then `msolve -P 1 -f flat8_Q.ms -o rur_Q.out`, `python3 verify_rur_exact.py`,
`python3 real_outside.py`, `Singular -q std_Q.sing`.
