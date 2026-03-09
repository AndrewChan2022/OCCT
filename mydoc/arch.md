# Architecture — OCCT Fitting & Approximation Stack

## Two Entry Paths

The fitting system has two fundamentally different entry paths depending
on whether the input is a continuous function or discrete points.

```
                        +-----------------------+
                        |     User / Caller     |
                        +-----------+-----------+
                                    |
              +---------------------+---------------------+
              |                                           |
     continuous function F(t)                  discrete points {P_i}
     (can evaluate at any t)                   (measured / sampled)
              |                                           |
   +----------v-----------+                  +------------v-----------+
   |   Approx  (wrappers) |                  |  AppDef  (application) |
   |   Approx_Curve3d     |                  |  AppDef_BSplineCompute |
   |   Approx_Curve2d     |                  |  AppDef_Compute        |
   |   Approx_CurveOnSurf |                  |  AppDef_Variational    |
   |   Approx_SweepApprox |                  +------------+-----------+
   +----------+-----------+                               |
              |                                           |
     +--------+--------+                        +---------v----------+
     |                  |                        |   AppParCurves     |
     v                  v                        |   (generic base)   |
 +---+------+    +------+------+                 |  LeastSquare.gxx   |
 | AdvApprox|    | AppCont     |                 |  Gradient.gxx      |
 | (1D core)|    | (simple LS) |                 |  BFGS optimizer    |
 +---+------+    +------+------+                 +---------+----------+
     |                  |                                  |
     v                  v                                  v
 +---+------+    +------+------+                 +---------+----------+
 |SimpleAppr|    | LeastSquare |                 |   MultiBSpCurve    |
 |Hermite+  |    | Gauss-      |                 |   (output)         |
 |Jacobi    |    | Legendre    |                 +--------------------+
 +---+------+    +-------------+
     |
     v
 +---+----------+
 | AdvApp2Var   |
 | (2D surface) |
 | tensor-prod  |
 +---+----------+
     |
     v
 +---+----------+
 | BSpline surf |
 | (output)     |
 +--------------+
```

## Module Dependency Graph

```
+-------------+     +---------------+     +------------------+
|  TKTopAlgo  |     |  TKGeomBase   |     |     TKG3d        |
|             |     |               |     |                  |
| BRepApprox -+---->| Approx -------+---->| AdvApprox        |
|             |     | AppDef -------+---->| (SimpleApprox,   |
|             |     | AppParCurves  |     |  ApproxAFunction)|
|             |     | AppCont       |     |                  |
|             |     | AdvApp2Var ---+---->|                  |
+-------------+     +---------------+     +------------------+
                           |                      |
                           v                      v
                    +------+------+        +------+------+
                    | math_BFGS   |        | PLib        |
                    | math_Matrix |        | Jacobi poly |
                    | FEmTool     |        | Hermite     |
                    +-------------+        | Gauss quad  |
                                           +-------------+
```

## Layered Architecture (bottom-up)

```
Layer 0 — Math foundations
+------------------------------------------------------------------+
|  PLib_JacobiPolynomial   PLib::HermiteInterpolate                |
|  Gauss-Jacobi quadrature tables   Convert_CompPolynomialToPoles  |
|  math_BFGS   math_Matrix   FEmTool (finite elements)            |
+------------------------------------------------------------------+

Layer 1 — Core fitting engines
+-------------------------------+  +-------------------------------+
|  AdvApprox_SimpleApprox       |  |  AppParCurves generics        |
|  (Hermite + Jacobi, 1 seg)   |  |  LeastSquare, Gradient, BFGS  |
+-------------------------------+  +-------------------------------+

Layer 2 — Multi-segment drivers
+-------------------------------+  +-------------------------------+
|  AdvApprox_ApproxAFunction    |  |  AppDef_BSplineCompute        |
|  (adaptive subdivision)       |  |  AppDef_Compute               |
|  AdvApp2Var_ApproxAFunc2Var   |  |  AppDef_Variational           |
|  (2D tensor-product patches)  |  |  AppCont_LeastSquare          |
+-------------------------------+  +-------------------------------+

Layer 3 — API wrappers
+-------------------------------+  +-------------------------------+
|  Approx_Curve3d               |  |  BRepApprox_Approx            |
|  Approx_Curve2d               |  |  (intersection curve fitting) |
|  Approx_CurveOnSurface        |  +-------------------------------+
|  Approx_SweepApproximation    |
+-------------------------------+

Layer 4 — High-level geometry API
+------------------------------------------------------------------+
|  GeomAPI_Interpolate          GeomAPI_PointsToBSpline            |
|  GeomAPI_PointsToBSplineSurface                                  |
+------------------------------------------------------------------+
```

## The Two Pipelines — Detailed Flow

### Pipeline A: Continuous Rebuild

```
Adaptor3d_Curve / Adaptor2d_Curve2d / SweepFunction
    |
    |  .D0(t)  .D1(t)  .D2(t)
    v
Approx_Curve3d / SweepApproximation
    |
    |  wraps as AdvApprox_EvaluatorFunction
    v
AdvApprox_ApproxAFunction
    |
    |  for each interval [a,b]:
    |      SimpleApprox.Perform([a,b])
    |          1. Hermite R(t) at endpoints
    |          2. Gauss-Jacobi -> Jacobi coefficients c_k
    |          3. P(t) = R(t) + W(t)*Q(t)
    |          4. Degree reduction
    |      if error > tol: cut and retry
    |
    |  PrepareConvert: verify joint continuity
    |  Convert_CompPolynomialToPoles: -> BSpline
    v
Geom_BSplineCurve / Geom_BSplineSurface
```

### Pipeline B: Discrete Point Fitting

```
Discrete points {P_i} with optional tangent/curvature
    |
    v
AppDef_MultiLine  (stores points as MultiPointConstraint rows)
    |
    |  parametrize: chord-length / centripetal / uniform
    v
AppDef_BSplineCompute
    |
    |  BFGS loop:
    |      1. Fix {u_i}, solve LS for poles Q
    |         (A^T A) Q = A^T P
    |      2. E = sum || C(u_i) - P_i ||^2
    |      3. dE/du_i = 2 (C(u_i)-P_i) . C'(u_i)
    |      4. BFGS update {u_i}
    |      5. Repeat
    |
    |  if cutting: split at max-error point, recurse
    v
AppParCurves_MultiBSpCurve
    |
    |  Approx_MCurvesToBSpCurve (join segments)
    v
Geom_BSplineCurve
```

### Pipeline C: Variational (Fair Curve)

```
Discrete points {P_i} + smoothness weights
    |
    v
AppDef_Variational
    |
    |  Minimize:
    |    E = w1*E_data + w2*||C'||^2 + w3*||C''||^2 + w4*||C'''||^2
    |
    |  Assemble FEM stiffness matrix
    |  Solve linear system for poles
    |  Iterate w1 if needed
    v
Geom_BSplineCurve  (smooth / fair)
```

## API Entry Points Summary

| Goal                                  | API Entry Point                          | Pipeline |
|---------------------------------------|------------------------------------------|----------|
| Rebuild analytic curve as BSpline     | `Approx_Curve3d`                         | A        |
| Rebuild curve-on-surface              | `Approx_CurveOnSurface`                 | A        |
| Rebuild sweep surface                 | `Approx_SweepApproximation`             | A        |
| Approximate 1D function -> BSpline    | `AdvApprox_ApproxAFunction`             | A        |
| Approximate 2D function -> surface    | `AdvApp2Var_ApproxAFunc2Var`            | A        |
| Fit discrete points -> BSpline curve  | `AppDef_BSplineCompute`                 | B        |
| Fit discrete points -> Bezier curve   | `AppDef_Compute`                        | B        |
| Fair/smooth curve through points      | `AppDef_Variational`                    | C        |
| Interpolate points (exact pass-thru)  | `GeomAPI_Interpolate`                   | (direct) |
| Interpolate surface grid              | `GeomAPI_PointsToBSplineSurface`        | (direct) |
| Fit BRep intersection curve           | `BRepApprox_Approx`                     | B        |

## Module Roles — One-Line Summary

| Module         | One-line role                                                |
|----------------|--------------------------------------------------------------|
| `AdvApprox`    | 1D Hermite+Jacobi engine with adaptive subdivision           |
| `AdvApp2Var`   | 2D tensor-product surface via repeated 1D AdvApprox          |
| `AppParCurves` | Generic LS + BFGS templates for parallel curve fitting       |
| `AppCont`      | Lightweight continuous-function Bezier LS (no BFGS)          |
| `AppDef`       | Application-level discrete-point fitter (uses AppParCurves)  |
| `Approx`       | User-facing wrappers that hide AdvApprox/AppCont behind APIs |
| `BRepApprox`   | BRep intersection fitting (AppParCurves with BRep types)     |
| `PLib`         | Jacobi polynomials, Hermite interpolation, Gauss quadrature  |
