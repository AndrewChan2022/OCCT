# AppDef — High-Level Multiline Approximation

## Architecture

```
AppDef_MultiLine                    (input: discrete multi-curve point data)
    |
    +-- AppDef_MultiPointConstraint  (one row: points + optional tangent + curvature)
    |       (inherits AppParCurves_MultiPoint)
    |
    +-- AppDef_MyLineTool            (adaptor: exposes MultiLine to generic templates)

Fitters (instantiated from AppParCurves generics):
    AppDef_Compute              (Bezier approximation, gradient/BFGS)
    AppDef_BSplineCompute       (BSpline approximation, gradient/BFGS)
    AppDef_Variational          (variational / fair-curve with regularization)
```

Internal generated classes (from .gxx templates):
```
AppDef_MyGradientOfCompute
AppDef_MyBSplGradientOfBSplineCompute
AppDef_MyGradientbisOfBSplineCompute
AppDef_ParLeastSquare*
AppDef_ParFunction*
AppDef_BSpParLeastSquare*
AppDef_BSpParFunction*
AppDef_Gradient_BFGS*
AppDef_BSpGradient_BFGS*
AppDef_LinearCriteria / AppDef_SmoothCriterion
```

## Purpose

`AppDef` is the **application-level approximation package**. It is the standard interface for fitting BSpline / Bezier curves to **discrete point sets with optional derivative constraints**.

Use cases:
- Fitting curves through measured data points (positions only, or with tangents/curvatures).
- Rebuilding a curve from sampled geometry.
- Fair-curve construction by variational energy minimization.

`AppDef` sits on top of `AppParCurves` (generic algorithms) and `AdvApprox` (adaptive subdivision); it wires them together and exposes friendly constructors.

## Principle / Theory / Algorithm

### Data Model

**`AppDef_MultiPointConstraint`** — one constraint at parameter `u`:
- 3D and/or 2D point positions.
- Optional tangent vectors (C1 constraint).
- Optional curvature vectors (C2 constraint, encoded as `κ·n` where `κ` is the curvature magnitude and `n` is the principal normal).

**`AppDef_MultiLine`** — ordered sequence of `MultiPointConstraint` objects indexed by parameter. This represents a set of parallel curves sampled at the same `u` values.

### Approximation Strategies

#### 1. `AppDef_Compute` — Bezier with BFGS parameter optimization

1. Initialize parameters `{uᵢ}` by chord-length / centripetal / uniform.
2. Solve least-squares for Bezier poles given current `{uᵢ}` (`AppDef_ParLeastSquare`).
3. Evaluate objective `E = Σ dist²`.
4. Compute gradient `∂E/∂uᵢ` (`AppDef_MyGradientOfCompute`).
5. Update `{uᵢ}` via BFGS (`AppDef_Gradient_BFGS*`).
6. Repeat until `E < tol²` or max iterations.
7. If tolerance not reached and `cutting = true`, subdivide at the point of maximum error and recurse.

#### 2. `AppDef_BSplineCompute` — BSpline with BFGS parameter optimization

Same loop as above but:
- Uses BSpline basis (`AppDef_BSpParLeastSquare`).
- Two sub-variants of the gradient:
  - `MyBSplGradient`: full BFGS on parameters with BSpline basis.
  - `MyGradientbis`: hybrid — uses Bezier locally then converts.
- Supports user-specified knots and multiplicities.
- `Interpol()` constructs an exact C2 cubic BSpline interpolation in one shot (no iterations).

#### 3. `AppDef_Variational` — Variational / Fair Curve

Minimises a combined objective:

```
E = λ₁ · E_data + λ₂ · E_length + λ₃ · E_curvature + λ₄ · E_torsion
```

where:
- `E_data` = sum of squared distances from data points to curve (fidelity).
- `E_length`, `E_curvature`, `E_torsion` = integrated squared first/second/third derivatives (fairness / smoothness).

The weights `λᵢ` are set via `AppDef_SmoothCriterion` / `AppDef_LinearCriteria`.

This is a linear problem in the BSpline poles (after fixing the knot vector) solved via `FEmTool` (finite-element assembly):
1. Assemble the stiffness matrix for each energy term over Hermite elements.
2. Add the data-fidelity matrix.
3. Solve the linear system for the optimal poles.
4. If the fidelity error is too large, increase `λ₁` and re-solve.

### Parametrization Types (`Approx_ParametrizationType`)

| Type | Formula |
|---|---|
| `ChordLength` | `uᵢ = uᵢ₋₁ + ‖Pᵢ - Pᵢ₋₁‖` |
| `Centripetal` | `uᵢ = uᵢ₋₁ + ‖Pᵢ - Pᵢ₋₁‖^0.5` |
| `IsoParametric` | `uᵢ = i / (n-1)` (uniform) |

## Files

### `AppDef_MultiPointConstraint.hxx / .cxx`
Data for one parameter value. Inherits `AppParCurves_MultiPoint`. Adds:
- `tabTang` / `tabTang2d`: tangent vectors.
- `tabCurv` / `tabCurv2d`: curvature vectors.
- `IsTangencyPoint()`, `IsCurvaturePoint()`: query flags.

### `AppDef_MultiLine.hxx / .cxx`
Container for `N` `MultiPointConstraint` objects. Supports:
- Construction from raw point arrays (no tangents).
- `SetParameter(Index, U)`: set the parameter value for row `Index`.
- `SetValue(Index, MPoint)`: set constraint data.

### `AppDef_MyLineTool.hxx / .cxx`
**Adaptor** that provides the interface expected by the `AppParCurves` generic templates:
- `FirstPoint()`, `LastPoint()`, `NbPoints()`.
- `Value(MLine, Index)` → `MultiPointConstraint`.
- `IsTangencyPoint(MLine, Index)`, `IsTheSamePoint(...)`.
- Used as a template parameter (not a base class).

### `AppDef_Compute.hxx / _0.cxx`
Bezier fitter. Key constructor parameters:
- `degreemin`, `degreemax`: degree search range.
- `Tolerance3d`, `Tolerance2d`.
- `NbIterations`: BFGS iterations.
- `cutting`: enable segment subdivision.
- `parametrization`: chord / centripetal / uniform.
- `Squares`: if true, skip BFGS and do direct least-squares only.

### `AppDef_BSplineCompute.hxx / _0.cxx`
BSpline fitter. Same interface as `Compute` plus:
- `SetKnots()`, `SetKnotsAndMultiplicities()`: fix knot vector.
- `SetContinuity(C)`: continuity of the spline (0, 1, or 2).
- `SetPeriodic(bool)`: periodic curve fitting.
- `Interpol(Line)`: exact cubic C2 interpolation.

### `AppDef_Variational.hxx / .cxx`
Variational fitter. Key parameters:
- `NbPoints`: data points.
- `MaxDegree`, `MaxSegment`: BSpline complexity limit.
- `Continuity` (`GeomAbs_Shape`).
- `WithMinMax`, `WithCutting`: additional options.
- `SmoothCriterion` weights (length, curvature, torsion).

### `AppDef_SmoothCriterion.hxx / .cxx` and `AppDef_LinearCriteria.hxx / .cxx`
Define the energy weights for the variational method. `LinearCriteria` allows separate weights for positions, tangents, and curvatures. `SmoothCriterion` is the base interface.

### Generated `*_0.cxx` files
Each `*_0.cxx` file is the C++ instantiation of one `.gxx` generic template with `AppDef` types. They contain no algorithmic code — just `#include` directives that trigger template expansion.

## Summary Table

| File | Role | Key Concept |
|---|---|---|
| `AppDef_MultiPointConstraint` | One-row input data | Position + tangent + curvature |
| `AppDef_MultiLine` | Full input dataset | Ordered rows of MultiPointConstraint |
| `AppDef_MyLineTool` | Generic adaptor | Template interface to MultiLine |
| `AppDef_Compute` | Bezier fitter | Least-squares + BFGS parameter opt. |
| `AppDef_BSplineCompute` | BSpline fitter | BSpline LS + BFGS, `Interpol()` |
| `AppDef_Variational` | Variational fitter | Energy minimisation, fair curves |
| `AppDef_SmoothCriterion` | Energy interface | Fairness weight definition |
| `AppDef_LinearCriteria` | Energy implementation | Position/tangent/curvature weights |
| `AppDef_*BFGS*` (generated) | BFGS optimizer instances | Parameter update loop |
| `AppDef_*LeastSquare*` (generated) | LS solver instances | Normal equations per segment |
| `AppDef_*Function*` (generated) | Objective function instances | E(params) evaluation |


## notes

1. call AdvApprox is for adaptive segment
2. call BFGS for parameter optimization
3. Variational add regularization term for smooth

## notes 2
They solve two completely different problems:

----
AppDef BFGS — optimizes the parameter u assigned to each input point

The input is discrete data points {P₁, P₂, …, Pₙ}. Parameters {u₁, u₂, …, uₙ} are not known — they must be estimated. The process is:


1. Init uᵢ by chord-length / centripetal / uniform
2. Fit BSpline poles Q given fixed {uᵢ}   ← least-squares (linear solve)
3. Compute E = Σ ||C(uᵢ) - Pᵢ||²
4. Compute ∂E/∂uᵢ                         ← gradient
5. Update {uᵢ} via BFGS                   ← nonlinear parameter optimization
6. Repeat 2–5 until convergence
BFGS is purely about finding better u values for the points. The poles are re-solved by least-squares at each BFGS step.

----
AdvApprox — adaptive segment subdivision for a continuous function

The input is an evaluable function F(t) (you can call it at any t). The problem is how many segments and where to cut:


1. Try to fit [a, b] with one polynomial (Hermite + Jacobi)
2. If error > tolerance → cut [a, b] → [a, c] + [c, b]
3. Repeat recursively until all segments are within tolerance

No parameter optimization — the parameters are just the quadrature points chosen internally.

----

side-by-side AppDef BFGS vs AdvApprox

| | AppDef BFGS | AdvApprox |
|---|---|---|
| Input | Discrete points `{Pᵢ}` | Continuous function `F(t)` |
| Unknown | Parameter `uᵢ` per point | Number + location of segments |
| Core math | Nonlinear optimization (BFGS) | Polynomial approximation + bisection |
| Sampling | Fixed (the N input points) | Adaptive (Gauss points per segment) |
| Continuity | Enforced by knot multiplicity choice | Enforced by Hermite endpoint constraints |
| Used by | AppDef, BRepApprox | Approx_Curve3d, Approx_SweepApproximation |

**AppDef BFGS** optimizes the parameter `uᵢ` assigned to each discrete input point:
1. Init `uᵢ` by chord-length / centripetal / uniform
2. Fit BSpline poles `Q` given fixed `{uᵢ}` ← least-squares
3. Compute `E = Σ ||C(uᵢ) - Pᵢ||²`
4. Compute `∂E/∂uᵢ` ← gradient
5. Update `{uᵢ}` via BFGS ← nonlinear parameter optimization
6. Repeat until convergence

**AdvApprox** adaptively subdivides for a continuous evaluable function `F(t)`:
1. Try to fit `[a, b]` with one polynomial (Hermite + Jacobi)
2. If error > tolerance → cut into `[a, c]` + `[c, b]`
3. Repeat until all segments are within tolerance

They are independent — `AppDef_BSplineCompute` does **not** use `AdvApprox` internally.