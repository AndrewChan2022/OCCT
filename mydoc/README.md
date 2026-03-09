# OCCT Fitting / Approximation — Documentation Index

## Module Map

```
TKG3d
└── AdvApprox          →  AdvApprox.md          1D adaptive function approximation (core engine)

TKGeomBase
├── AdvApp2Var         →  AdvApp2Var.md          2D surface approximation (extends AdvApprox)
├── AppParCurves       →  AppParCurves.md        Generic base: multi-curve least-squares + BFGS
├── AppCont            →  AppCont.md             Continuous-function Bezier fit (single segment)
├── AppDef             →  AppDef.md              Application-level discrete-point fitter
└── Approx             →  Approx.md              User-facing wrappers (Curve3d, Sweep, etc.)

TKTopAlgo
└── BRepApprox         →  TKTopAlgo_BRepApprox.md  BRep intersection curve fitting
```

## Layered Architecture

```
User code / geometry algorithms
        │
   Approx_Curve3d / Approx_SweepApproximation
        │                         │
   AdvApprox_ApproxAFunction   AdvApp2Var_ApproxAFunc2Var
        │                         │
   AdvApprox_SimpleApprox      (calls AdvApprox Iso + Patch logic)
        │
   PLib_JacobiPolynomial  (orthogonal polynomial basis)
        │
   Gauss–Jacobi quadrature  (integration → coefficients)

AppDef_BSplineCompute  ─────┐
AppDef_Compute         ─────┤── AppParCurves (generic templates)
AppDef_Variational     ─────┘     ├── LeastSquare (normal equations)
                                  ├── Function (objective E)
                                  ├── Gradient (∂E/∂u)
                                  └── BFGS (parameter update)
BRepApprox_Approx  ─── same AppParCurves templates, BRep MultiLine types
```

## Demo Files

Located in `mydemo/`:

| File | Demos |
|---|---|
| [demo_interp.cpp](../mydemo/demo_interp.cpp) | 1–4: Interpolation (points, end tangents, per-point tangents, surface grid) |
| [demo_approx.cpp](../mydemo/demo_approx.cpp) | 5–8: Approximation rebuild (single curve, chain, discontinuity, surface) |

## Key Entry Points for New Code

| Goal | API |
|---|---|
| Fit curve through measured 3D points | `GeomAPI_Interpolate` |
| Fit curve through 3D points with tangents | `GeomAPI_Interpolate::Load(tangents, flags)` |
| Fit surface through grid points | `GeomAPI_PointsToBSplineSurface::Interpolate()` |
| Rebuild any analytic curve as BSpline | `Approx_Curve3d(adaptor, tol, C2, maxSeg, maxDeg)` |
| Rebuild sweep surface | `Approx_SweepApproximation::Perform()` |
| Fit discrete points → BSpline (with BFGS) | `AppDef_BSplineCompute` |
| Fair curve (variational energy) | `AppDef_Variational` |
| Fit 2D function → BSpline surface | `AdvApp2Var_ApproxAFunc2Var` |
| Fit any 1D function → BSpline | `AdvApprox_ApproxAFunction` |
