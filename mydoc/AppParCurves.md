# AppParCurves — Parallel Curve Approximation Base Layer

## Architecture

```
AppParCurves               (utility: Bernstein / BSpline matrix helpers)
    |
    +-- AppParCurves_MultiPoint        (one column: N points at the same parameter)
    |
    +-- AppParCurves_MultiCurve        (result: N parallel Bezier curves)
    |
    +-- AppParCurves_MultiBSpCurve     (result: N parallel BSpline curves)
    |
    +-- AppParCurves_Constraint        (enum: PassPoint / TangencyPoint / CurvaturePoint)
    +-- AppParCurves_ConstraintCouple  (pair: index + constraint type)
    |
    +-- Generic templates (.gxx files, instantiated by AppDef / BRepApprox)
         AppParCurves_LeastSquare.gxx
         AppParCurves_Function.gxx
         AppParCurves_Gradient.gxx
         AppParCurves_BSpFunction.gxx
         AppParCurves_BSpGradient.gxx
         AppParCurves_ResolConstraint.gxx
```

## Purpose

`AppParCurves` is the **generic algorithmic base layer** for approximating a set of parallel discrete point sequences (a "MultiLine") by a set of parallel curves (a "MultiCurve").

It does not operate on real geometry directly. Instead, it provides:
- **Data structures**: `MultiPoint`, `MultiCurve`, `MultiBSpCurve`.
- **Numerical algorithms**: least-squares matrix assembly, gradient computation, BFGS optimization, constraint resolution — all as C++ generic templates (`.gxx` files).
- **Utility functions**: Bernstein matrix, BSpline matrix, second derivatives.

Concrete packages (`AppDef`, `BRepApprox`) instantiate these generics with their own `MultiLine` and `LineTool` types.

## Principle / Theory / Algorithm

### 1. The MultiLine Model

Input is conceptually a matrix of points:

```
Parameter   u₁     u₂     ...    uₙ
Curve 1:    P₁¹    P₁²    ...    P₁ⁿ
Curve 2:    P₂¹    P₂²    ...    P₂ⁿ
  ...
Curve k:    Pₖ¹    Pₖ²    ...    Pₖⁿ
```

All `k` curves share the same parameter sequence `{uᵢ}`.

### 2. Least-Squares Fitting

Given parameters `{uᵢ}`, fit `k` Bezier (or BSpline) curves simultaneously by solving one least-squares system. Because each curve is independent in its coordinates but shares the same basis matrix `A(u)`, the system decouples per coordinate component:

```
min Σᵢ ||C(uᵢ) - Pᵢ||²
```

The normal equations become:

```
(AᵀA) Q = AᵀP
```

where `A` is the Bernstein (or BSpline) collocation matrix and `Q` are the poles. All curves share the same `AᵀA` matrix, only the right-hand side `AᵀP` differs per curve / per dimension.

### 3. Parametrization

Parameters `{uᵢ}` are not given directly by the caller — they are computed (chord-length, centripetal, or uniform) and **iteratively refined** by the gradient / BFGS optimizer to minimise the total residual.

### 4. Gradient / BFGS Optimization

`AppParCurves_Function.gxx`: defines the objective function `E(u₁,…,uₙ)` (sum of squared distances from each point to its nearest point on the curve at the current parameter).

`AppParCurves_Gradient.gxx`: computes `∂E/∂uᵢ` analytically using the curve derivative.

`AppParCurves_BSpGradient.gxx` / `BSpFunction.gxx`: BSpline variants.

The BFGS optimizer (`math_BFGS` wrapped by the `_BFGS` classes) iterates: recompute parameters → refit → recompute gradient → update until convergence.

### 5. Constraint Resolution

`AppParCurves_ResolConstraint.gxx`: when endpoint tangency or curvature constraints are imposed, the constrained poles are solved via a bordered system. The free poles are computed from the least-squares normal equations after removing the constrained poles' contributions.

## Files

### `AppParCurves.hxx / .cxx`
**Utility class.**
- `BernsteinMatrix(NbPoles, U, A)`: fills collocation matrix for Bezier evaluation.
- `Bernstein(NbPoles, U, A, DA)`: fills value and first-derivative Bernstein matrices.
- `SecondDerivativeBernstein(U, DDA)`: second derivatives (for curvature constraints).
- `SplineFunction(NbPoles, Degree, Parameters, FlatKnots, A, DA, Index)`: BSpline collocation matrix.

### `AppParCurves_MultiPoint.hxx / .cxx / .lxx`
One column of the MultiLine. Stores a mixed set of 3D and 2D points at the same parameter value. No constraint data — that lives in `AppDef_MultiPointConstraint`.

### `AppParCurves_MultiCurve.hxx / .cxx`
Result of a Bezier approximation: a set of parallel Bezier curves encoded as their control poles. Provides `Value()`, `D1()`, `D2()` for evaluation.

### `AppParCurves_MultiBSpCurve.hxx / .cxx`
Result of a BSpline approximation: parallel BSpline curves with a shared knot vector and multiplicities.

### `AppParCurves_Constraint.hxx`
Enum: `NoConstraint`, `PassPoint`, `TangencyPoint`, `CurvaturePoint`.

### `AppParCurves_ConstraintCouple.hxx / .cxx`
Pair `(pointIndex, Constraint)`. Used to record which data points have constraints.

### `AppParCurves_LeastSquare.gxx`
Generic least-squares solver for one segment. Assembles the Bernstein collocation matrix and solves the normal equations with LAPACK (via `math_Matrix`). Handles endpoint constraints by partitioning the system.

### `AppParCurves_Function.gxx`
Generic objective function (Bezier variant). Evaluates `E` and provides the parameter-update step for gradient descent.

### `AppParCurves_Gradient.gxx`
Gradient of `E` with respect to the parameters. Used by the BFGS loop.

### `AppParCurves_BSpFunction.gxx` / `BSpGradient.gxx`
BSpline variants of `Function` and `Gradient`.

### `AppParCurves_ResolConstraint.gxx`
Resolves endpoint constraints (tangency, curvature) by partitioning the pole vector into free and constrained parts and solving the reduced system.

## Summary Table

| File | Role | Key Concept |
|---|---|---|
| `AppParCurves.hxx/.cxx` | Matrix utilities | Bernstein / BSpline collocation |
| `AppParCurves_MultiPoint` | Column of input data | Mixed 2D/3D points at same `u` |
| `AppParCurves_MultiCurve` | Bezier result | Parallel Bezier control poles |
| `AppParCurves_MultiBSpCurve` | BSpline result | Parallel BSpline with shared knots |
| `AppParCurves_Constraint.hxx` | Constraint enum | Pass / Tangency / Curvature |
| `AppParCurves_ConstraintCouple` | Index + type pair | Per-point constraint record |
| `AppParCurves_LeastSquare.gxx` | Bezier least-squares | Normal equations, constrained solve |
| `AppParCurves_Function.gxx` | Bezier objective | E(params), parameter update |
| `AppParCurves_Gradient.gxx` | Bezier gradient | ∂E/∂u analytic gradient |
| `AppParCurves_BSpFunction.gxx` | BSpline objective | Same as Function but BSpline |
| `AppParCurves_BSpGradient.gxx` | BSpline gradient | ∂E/∂u for BSpline |
| `AppParCurves_ResolConstraint.gxx` | Constraint resolution | Bordered system for constrained poles |
