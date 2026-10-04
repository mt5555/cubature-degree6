# Degree 9: numerical search only (no proof)

**Status: open.** The best known rule with all nodes inside has 19 points (fully symmetric;
Lyness–Jespersen type, listed in Cools' encyclopedia and by Wandzura–Xiao). Papanicolopulos
(2015) proved 19 is minimal among *fully symmetric* rules. The lower bound is 17: the basic
bound dim P₄ = 15, raised by a Möller-type argument, since Gaussian rules don't exist on the
triangle. For rules without full symmetry, 17 and 18 points were open as far as we know.

These are **numerical findings only**. They are evidence, not a proof.

## Method

`search9.py N mode starts seed` runs node-space least squares with variable projection: the
weights are the linear least-squares solution given the nodes. The residual is measured in
orthonormal-polynomial coordinates over all 55 monomials of degree ≤ 9, where the zero rule
has residual 0.707. In mode `inside`, nodes are parametrized by squared barycentric
coordinates. In mode `free`, nodes are unconstrained.

## Results

| N | mode | starts | exact rules (residual < 1e-10) | best residual |
|---|---|---|---|---|
| 19 | inside | 30 | 15 (13 with positive weights, all inside) | 9×10⁻¹² |
| 18 | inside | 50 | 0 | 0.0189 (same local minimum from every start) |
| 18 | free | 50 | 0 | 0.0050 |
| 17 | inside | 150 | 0 | 0.138 |
| 17 | free | 50 | 0 | 0.108 |

* The N = 19 runs check that the method works: half the random starts give an exact rule.
* N = 17 and 18 are overdetermined (51 and 54 unknowns against 55 equations), so a solution would
  need special structure. None was found, even with nodes allowed outside.
* The best 18-point free attempts send 2–4 nodes 5–8 units outside with weights about 10⁻¹⁴.
  So they effectively use fewer nodes, and still don't reach an exact rule.
* A proof would be much harder than at degrees 6 and 8. For N > dim P₄ the moment matrix is
  not a flat extension at the first level, so there are many more unknowns.

Result files: `r9_N*_*.npy`, with one row per start: [residual, max edge violation, min weight, x, y, w].
