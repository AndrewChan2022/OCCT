# AppCont — Continuous Function Approximation

## Architecture

```
AppCont_Function          (abstract: user-supplied continuous function)
    |
    +-- AppCont_LeastSquare   (single-segment Bezier least-squares fitter)
```

`AppCont_ContMatrices*.pxx` — pre-computed polynomial matrices embedded as data.

## Purpose

`AppCont` approximates a **continuous function** (as opposed to `AppDef`'s discrete point sets) by a set of parallel Bezier curves in one shot using least-squares.

The user provides a concrete sub-class of `AppCont_Function` with methods to evaluate the function and its first derivative at any parameter value. `AppCont_LeastSquare` samples that function at Gauss–Legendre points and solves the normal equations directly without iterative parameter refinement.

This module is simpler than `AppDef`/`AppParCurves` because the parameters are fixed (Gauss points — not optimised), so only a single linear solve is needed.

## Principle / Theory / Algorithm

### 1. Gauss Sampling

`AppCont_LeastSquare` selects `NbPoints` Gauss–Legendre quadrature abscissae on `[U0, U1]`. The function is evaluated at these points (via `AppCont_Function::Value()` and `D1()`).

### 2. Collocation Matrix

The Bezier basis matrix `B` of degree `Deg` is assembled at the Gauss points using the utility routines from `AppParCurves` (Bernstein matrix).

### 3. Least-Squares Solve

The normal equations

```
BᵀB · Q = Bᵀ · P
```

are solved for the poles `Q`. Because the Gauss points are fixed, `BᵀB` depends only on degree and can in principle be pre-computed (hence the `ContMatrices*.pxx` pre-computed tables).

### 4. Endpoint Constraints

If `FirstCons` or `LastCons` specifies `TangencyPoint` or `CurvaturePoint`, the first/last 1–2 poles are fixed from the function's derivative values at the endpoints, and the system is solved for the remaining free poles.

### 5. Periodicity Support

`PeriodicityInfo` in `AppCont_LeastSquare` allows one output dimension to be treated as an angular parameter modulo `2π` (e.g., for curves on surfaces of revolution).

## Files

### `AppCont_Function.hxx`
**Abstract base class.**
Users must implement:
- `FirstParameter()` / `LastParameter()` — domain.
- `Value(u, Pnt2d[], Pnt[])` — evaluate function values.
- `D1(u, Vec2d[], Vec[])` — evaluate first derivatives.
- `PeriodInformation(dimIdx, isPeriodic, period)` — optional periodicity hint.
- Members `myNbPnt` (3D outputs) and `myNbPnt2d` (2D outputs) must be set.

### `AppCont_LeastSquare.hxx / .cxx`
**Least-squares solver for a single Bezier segment.**
- Constructor: takes the `AppCont_Function`, parameter range, endpoint constraints, degree, and number of sample points.
- `Value()`: returns the fitted `AppParCurves_MultiCurve` (parallel Bezier curves).
- `Error(F, MaxE3d, MaxE2d)`: returns RMS error `F`, max 3D error, max 2D error.
- `IsDone()`: whether the solve succeeded.

Internal layout:
- `myPoints`: matrix of sampled function values (rows = sample points, cols = total dimension).
- `myVB`: Bernstein matrix at sample parameters.
- `myPoles`: solved control poles.

### `AppCont_ContMatrices*.pxx`
Pre-computed data for the `BᵀB` normal equation matrix at standard degrees and Gauss point counts:
- `ContMatrices_BB.pxx`: `BᵀB` matrices.
- `ContMatrices_IBP.pxx`: `(BᵀB)⁻¹ Bᵀ` projection matrices (for unconstrained solve).
- `ContMatrices_IBT.pxx`: projected matrices with tangency constraints.
- `ContMatrices_InvM.pxx`: inverse matrices.
- `ContMatrices_VB.pxx`: Bernstein values at Gauss points.

## Summary Table

| File | Role | Key Concept |
|---|---|---|
| `AppCont_Function.hxx` | Abstract evaluator | Continuous function oracle D0+D1 |
| `AppCont_LeastSquare` | Bezier fit (single segment) | Gauss sampling + normal equations |
| `AppCont_ContMatrices*.pxx` | Pre-computed matrices | `BᵀB`, `(BᵀB)⁻¹Bᵀ` at standard degrees |
