# Approx Module — Rebuild (AdvApprox) vs Point Fitting (AppDef)

## Uses AdvApprox (continuous function → adaptive segment rebuild)

| Class | What it rebuilds |
|---|---|
| `Approx_Curve3d` | Any `Adaptor3d_Curve` → `Geom_BSplineCurve` |
| `Approx_Curve2d` | Any `Adaptor2d_Curve2d` → `Geom2d_BSplineCurve` |
| `Approx_CurveOnSurface` | Curve on surface → 3D + 2D PCurve simultaneously |
| `Approx_SweepApproximation` | Section-law sweep → BSpline surface |
| `Approx_SameParameter` | Reparametrization map between two curves |
| `Approx_CurvilinearParameter` | Arc-length reparametrization function |

All of these wrap a **continuous evaluable function** and call `AdvApprox_ApproxAFunction` directly.

## Uses AppCont (continuous function, direct LS per segment, no BFGS)

| Class | What it does |
|---|---|
| `Approx_FitAndDivide` | Fit + subdivide loop, uses `AppCont_Function` as the MultiLine type |
| `Approx_FitAndDivide2d` | Same, 2D variant |

These use the `AppParCurves_LeastSquare` machinery (like AppDef) but with `AppCont_Function` as
input — a **continuous function** rather than discrete points. No BFGS — just direct least-squares
at Gauss points per segment, then subdivide if error is too large.

## Nothing in `Approx` calls `AppDef` directly

`AppDef` is used upstream (by the caller) when the input is discrete points. Once you have a
continuous adaptor, you go through `AdvApprox` or `AppCont`.

## Summary Diagram

```
Input: discrete points {Pᵢ}
    └─→ AppDef_BSplineCompute  (BFGS parameter opt.)
         └─→ AppParCurves generics

Input: continuous function / adaptor
    └─→ Approx_Curve3d / Curve2d / CurveOnSurface
         └─→ AdvApprox_ApproxAFunction  (Hermite+Jacobi, adaptive segments)

Input: continuous function, Bezier segment per piece
    └─→ Approx_FitAndDivide / FitAndDivide2d
         └─→ AppCont_LeastSquare  (Gauss sample, direct LS, subdivide)
```
