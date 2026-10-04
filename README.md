# No 10-point degree-6 cubature rule on the triangle has all its nodes inside

**Theorem.** Let T = {x ≥ 0, y ≥ 0, x + y ≤ 1}. No cubature rule with 10 nodes, all in the
closed triangle T, integrates every polynomial of total degree ≤ 6 exactly. This holds for
any real weights. Since 10 = dim P₃ is also the lower bound, the minimal number of nodes for
a degree-6 rule with all nodes inside T is 11, the count achieved by Day & Taylor (PAMM 7, 2007).
Section 5 certifies rigorously that an 11-point rule with all nodes strictly inside exists.

This is a computer-assisted proof. Every step of the final certificate is checked in exact
rational arithmetic. Because exactness forces the nodes to be distinct (Section 1), the result
also rules out rules with fewer than 10 distinct nodes.

## 1. Reduction to 8 unknowns

Let μ = Σ wᵢ δ(xᵢ) be a 10-point rule that is exact to degree 6.

* Write V for the 10×10 Vandermonde matrix on P₃. Then Vᵀ diag(w) V = M₃, the Lebesgue moment
  matrix on P₃, which is positive definite. So V is invertible, the nodes are distinct and
  P₃-unisolvent, and every weight wᵢ is positive.
* Let M₄(μ) be the moment matrix of μ on P₄. It has rank 10 = rank M₃, so it is a *flat
  extension* of M₃. The only moments it involves that are not already fixed by exactness are
  the 8 degree-7 moments yₖ = μ(x^(7−k) yᵏ), k = 0..7. Flatness forces the degree-8 block to
  be C = Bᵀ M₃⁻¹ B. That block must also have moment (Hankel) structure, which gives
  **6 quadratic equations h(y) = 0**.
* For each g in {x, y, 1−x−y}, the localizing matrix M₃(g·μ) = Σ wᵢ g(xᵢ) v(xᵢ)v(xᵢ)ᵀ is
  positive semidefinite when every node is in the closed triangle. A node on the boundary only
  makes it singular. Its entries are moments of degree ≤ 7. The P₂ block A_g and the
  off-diagonal block F_g involve only moments of degree ≤ 6, so they are fixed, and A_g ≻ 0
  (checked exactly). Taking the Schur complement on A_g gives
  **three 4×4 linear matrix inequalities G_g(y) = D_g(y) − F_gᵀA_g⁻¹F_g ⪰ 0.** The block D_g(y) is
  affine in y: Hank(y₀..y₆) for g = x, Hank(y₁..y₇) for g = y, and
  Hank(m₆) − Hank(y₀..y₆) − Hank(y₁..y₇) for g = 1−x−y.

So an inside rule implies a point y ∈ ℝ⁸ in K ∩ V, where K = {G_g ⪰ 0} and V = {h = 0}. The
converse also holds (Curto–Fialkow flat extension theorem), but only this direction is needed here.

## 2. Certificate that K ∩ V is empty

Coordinates: tₖ = (yₖ − cₖ)/hwₖ, where c and hw are exact dyadic rationals and hw ≠ 0. Each LMI
is congruence-scaled by a dyadic rational matrix W_g (whitening): G ⪰ 0 ⇒ W G Wᵀ ⪰ 0, and W_g
is in fact nonsingular. Each hₖ is scaled by a nonzero rational constant, giving h̃ₖ.

**(a) Bounding box.** For each k and each sign there are exact rational matrices Z_g ≻ 0 with
±tₖ + c₀ = Σ_g ⟨Z_g, G̃_g(t)⟩ identically in t. So |tₖ| ≤ 1.03 on K. The linear identities
are made exact by an exact correction step. Positive definiteness is checked by exact LDLᵀ.

**(b) Main identity** (degree 4 in t, 495 monomials):

    P(t) := σ₀(t) + Σ_g ⟨S_g, G̃_g(t) ⊗ v₁v₁ᵀ⟩ + Σₖ λₖ(t) h̃ₖ(t) + γ

* σ₀ = v₂ᵀ S₀ v₂, where v₂ is the 45 monomials of degree ≤ 2. S₀ is the numerical Gram
  matrix plus 10⁻⁸·I, and it is exactly positive definite.
* Each S_g is a 36×36 exactly positive-definite rational matrix, and v₁ = (1, t₁, …, t₈).
* Each λₖ is a rational quadratic polynomial, and γ = 6.844×10⁻³.

On K ∩ V every term except γ is ≥ 0 (a Kronecker product of PSD matrices paired with a PD
matrix is ≥ 0) and h̃ = 0, so P(t\*) ≥ γ. But in exact arithmetic,
Σ |P_μ| · Πᵢ Bᵢ ≤ **5.0×10⁻⁷** < γ on the certified box ⊇ K. That is a contradiction, so
K ∩ V = ∅. ∎ (Nearly all of the 5.0×10⁻⁷ comes from the δ-shift. The rounding residual of the
identity itself is about 10⁻¹⁰.)

The certificate has to be shifted by a tiny δ·I, instead of rounded with an exact projection,
because the top-degree parts of the 6 equations all vanish on the curve yₖ = zᵏ. That forces
a 15-dimensional null space onto any exact Gram matrix S₀.

## 3. Supporting numerics (not part of the proof)

* The search tools reproduce the 12-point Dunavant rule and find 11-point inside rules from
  200 of 200 random starts. 10-point rules with nodes outside are easy to find.
* With all 10 nodes inside, the best residual found is 0.0294 out of 0.707. This was confirmed with
  constrained SLSQP: one node sits at a vertex, one on an edge.
* The best edge offset found that allows an exact 10-point rule is ε ≈ 0.10454. This is the best found, not certified.
* After scaling so that Lebesgue measure gives every eigenvalue = 1, the best value of the
  worst eigenvalue over exact 10-point rules is ≈ −0.934 (it would need to be ≥ 0). This comes
  from sampling plus local ascent over about 1,200 rules. The level-2 moment relaxation stays
  infeasible down to −0.3, which is why the degree-4 certificate exists.

## 4. Audits and independent review

* An independent referee (a separate Claude Fable 5.1 agent) wrote its own code to
  re-derive the 10×10 localizing matrices and the flatness equations exactly. It confirmed that
  the Schur forms, the 6 equations, and the certificate data match exactly. It also
  re-implemented the final check and the 16 box checks. Verdict: the proof checks out.
* `verify_data.py` (adapted from the referee's script) repeats the exact data re-derivation.
* The 4×4 Schur form matches the direct 10×10 localizing matrices to 2×10⁻¹⁰ at random points.
* Negative control: at an exact outside rule the certificate terms sum to γ + O(10⁻⁸). The
  (1−x−y) term is the one that goes negative there, as it should.

## 5. Certified 11-point rule with all nodes inside (existence)

`cert11.py` certifies that a genuine 11-point degree-6 rule exists with every node strictly
inside T. It uses the Krawczyk interval-Newton test in exact rational interval arithmetic.

* There are 33 unknowns and 28 equations. Pivoted QR picks 5 unknowns to fix (x₃, x₄, x₆, y₃,
  y₁₁), set to exact dyadic rationals. Newton's method in mpmath (70 digits) solves for the other
  28, reaching residual 10⁻⁷².
* Box X = center ± 10⁻³⁰. With an approximate inverse Y (exact rationals), the test
  K(X) = x̃ − Y·F(x̃) + (I − Y·J(X))(X − x̃) ⊂ int X holds. The worst ratio |Kᵢ − x̃ᵢ| / r is
  7.6×10⁻¹¹. So X contains exactly one exact rule.
* Over all of X, every node satisfies x, y, 1−x−y > 0 with edge distance ≥ 0.0453, and every
  weight is > 0 (the smallest is 0.016, with the weights summing to 1/2).
* The rule was chosen to maximize the minimum distance from a node to an edge. It is
  mirror-symmetric across x = y to about 10⁻¹³. It is not necessarily Day & Taylor's rule.

| node | x | y | weight (sum = 1/2) |
|---|---|---|---|
| 1 | 0.04534843557757456 | 0.04534843557759772 | 0.01607025470980284 |
| 2 | 0.07492814583849801 | 0.61100718692760725 | 0.05890436798786112 |
| 3 | 0.17310639355862059 | 0.17310639355866236 | 0.03859898544711080 |
| 4 | 0.63454745425993975 | 0.29646524581443862 | 0.05528316171582048 |
| 5 | 0.30271617543506396 | 0.04534843557764037 | 0.03594773427991819 |
| 6 | 0.61100718692757461 | 0.07492814583848484 | 0.05890436798785817 |
| 7 | 0.88053899089008947 | 0.05532863648360303 | 0.02355269993193612 |
| 8 | 0.29646524581529349 | 0.63454745425907988 | 0.05528316171587405 |
| 9 | 0.33513321853409461 | 0.33513321853407335 | 0.09795483201186181 |
| 10 | 0.05532863648392812 | 0.88053899088976817 | 0.02355269993204220 |
| 11 | 0.04534843557763532 | 0.30271617543513590 | 0.03594773427991421 |

**Day & Taylor's rule, certified.** The 11-point rule published in Day & Taylor (PAMM 7, 2007)
is in `daytaylor11.txt`. In floating point it has residual 4×10⁻¹³, minimum edge distance 0.0233
and minimum weight 0.019. Running `python3 cert11.py daytaylor11.npy cert11_daytaylor.pkl`
certifies it the same way. The fixed unknowns are x₃, x₅, x₆, y₆ and y₇. The Krawczyk ratio is
8×10⁻¹¹, and the certified exact rule lies within 10⁻¹⁵ of the published values. So the
published rule is a genuine 11-point degree-6 rule with all nodes strictly inside, and with
the proof above it is optimal.

## 6. Reproducing

Files are in `tri/`. These checks are standalone. They need only Python with numpy and sympy
(mpmath and scipy for cert11.py), and no SDP solver:

* `python3 verify_data.py` (about 1 min): exact re-derivation of the problem data from the moment matrices.
* `python3 rigcert.py` (about 3 s): exact re-check of the 16 box certificates and the main certificate.
  It uses `sosdata.pkl`, `rigbox.pkl` and `sos_num.pkl`.
* `python3 cert11.py` (about 3 s): Krawczyk certification of the 11-point rule. It uses `best11_1.npy`.
  Add the arguments `daytaylor11.npy cert11_daytaylor.pkl` to certify Day & Taylor's published rule instead.

To rebuild from scratch, run `box.py` → `sosdata.py` → `sos.py` → `rigbox.py` → `rigcert.py`
(this needs cvxpy with Clarabel). The supporting and exploratory scripts are `common.py`,
`search.py`, `minout.py`, `refine.py`, `cf.py`, `cf2.py`, `fixA*.py`, `lasserre.py`,
`best11.py` and `audit.py`.
