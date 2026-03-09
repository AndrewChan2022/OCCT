# AppDef vs AdvApprox — Two Separate Fitting Stacks

**AppDef has its own fitting — it does NOT call AdvApprox.** They are completely separate pipelines.

## Two Completely Separate Stacks

```
Discrete points {Pᵢ}          Continuous function F(t)
       │                               │
  AppDef_BSplineCompute          Approx_Curve3d
       │                               │
  AppParCurves generics          AdvApprox_ApproxAFunction
       │                               │
  LeastSquare (normal eqns)      SimpleApprox (Hermite+Jacobi)
  Gradient (∂E/∂u)               DichoCutting / PrefAndRec
  BFGS (update u)                Convert_CompPolynomialToPoles
       │                               │
  MultiBSpCurve                  BSpline poles + knots
```

## AppDef BFGS — Optimizes the Parameter `u` Assigned to Each Input Point

The input is **discrete data points** `{P₁, P₂, …, Pₙ}`. Parameters `{u₁, u₂, …, uₙ}` are not
known — they must be estimated. The process is:

1. Init `uᵢ` by chord-length / centripetal / uniform
2. Fit BSpline poles `Q` given fixed `{uᵢ}` ← least-squares (linear solve)
3. Compute `E = Σ ||C(uᵢ) - Pᵢ||²`
4. Compute `∂E/∂uᵢ` ← gradient
5. Update `{uᵢ}` via BFGS ← nonlinear parameter optimization
6. Repeat 2–5 until convergence

BFGS is purely about finding better `u` values for the points. The poles are re-solved by
least-squares at each BFGS step.

## AdvApprox — Adaptive Segment Subdivision for a Continuous Function

The input is an **evaluable function** `F(t)` (you can call it at any `t`). The problem is
how many segments and where to cut:

1. Try to fit `[a, b]` with one polynomial (Hermite + Jacobi)
2. If error > tolerance → cut `[a, b]` → `[a, c]` + `[c, b]`
3. Repeat recursively until all segments are within tolerance

No parameter optimization — the parameters are just the quadrature points chosen internally.

## AppDef's Own "Subdivision"

`AppDef_BSplineCompute` does have a `cutting` option, but it is **not** the same as AdvApprox's
adaptive segmentation. It is just:

```
if error > tolerance AND cutting == true:
    split the point set at the worst point
    re-run the whole BFGS fit on each half separately
```

No Hermite, no Jacobi, no Gauss points. Just a re-fit on a smaller point subset.

## Side-by-Side Comparison

| | AppDef BFGS | AdvApprox |
|---|---|---|
| Input | Discrete points `{Pᵢ}` | Continuous function `F(t)` |
| Unknown | Parameter `uᵢ` per point | Number + location of segments |
| Core math | Nonlinear optimization (BFGS) | Polynomial approximation + bisection |
| Sampling | Fixed (the N input points) | Adaptive (Gauss points per segment) |
| Continuity | Enforced by knot multiplicity choice | Enforced by Hermite endpoint constraints |
| Used by | AppDef, BRepApprox | Approx_Curve3d, Approx_SweepApproximation |

They are independent — `AppDef_BSplineCompute` does **not** use `AdvApprox` internally.

## Which to Use When

| Situation | Use |
|---|---|
| You have measured/sampled discrete points | `AppDef_BSplineCompute` |
| You have a mathematical function you can evaluate anywhere | `Approx_Curve3d` → `AdvApprox` |
| You have a CAD curve and want a BSpline approximation | `Approx_Curve3d` (wraps the adaptor) |
| You want a fair/smooth curve through scattered points | `AppDef_Variational` |
