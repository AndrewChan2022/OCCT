# Approx — High-Level Approximation Wrappers

## Architecture

```
Approx (user-facing wrappers)
    |
    +-- Approx_Curve3d              wraps Adaptor3d_Curve → Geom_BSplineCurve
    +-- Approx_Curve2d              wraps Adaptor2d_Curve2d → Geom2d_BSplineCurve
    +-- Approx_CurveOnSurface       wraps a curve-on-surface → 3D + 2D BSpline
    +-- Approx_SweepApproximation   section-law sweep → BSpline surface
    |       uses AdvApprox_ApproxAFunction internally
    |
    +-- Approx_MCurvesToBSpCurve    merge multiple Bezier curves into one BSpline
    +-- Approx_SameParameter        reparametrize a curve to match another
    +-- Approx_CurvilinearParameter reparametrize by arc length
    +-- Approx_CurvlinFunc          arc-length function helper
    |
    +-- Generic templates (instantiated by other packages)
         Approx_ComputeLine.gxx        Bezier multi-curve fit
         Approx_ComputeCLine.gxx       continuous-function multi-curve fit
         Approx_BSplComputeLine.gxx    BSpline multi-curve fit
         Approx_FitAndDivide.gxx       fit + subdivide loop (3D)
         Approx_FitAndDivide2d.gxx     fit + subdivide loop (2D)
```

## Purpose

`Approx` is the **final user-facing layer** of the OCCT approximation stack. It converts geometric adaptor objects (`Adaptor3d_Curve`, etc.) or section-law functions into standard OCCT geometry types (`Geom_BSplineCurve`, `Geom_BSplineSurface`).

It hides the details of `AdvApprox`, `AppParCurves`, and `AppDef` behind simple constructors.

## Principle / Theory / Algorithm

### `Approx_Curve3d` and `Approx_Curve2d`

These call `AdvApprox_ApproxAFunction` with an internal evaluator that wraps the adaptor's `D0`/`D1`/`D2` methods. The algorithm is fully adaptive (see `AdvApprox` doc).

Steps:
1. Create an `EvaluatorFunction` that calls `Adaptor3d_Curve::D0/D1/D2`.
2. Call `AdvApprox_ApproxAFunction` with tolerance, continuity, degree and segment limits.
3. Convert the resulting poles and knots into a `Geom_BSplineCurve`.

### `Approx_CurveOnSurface`

Approximates a curve lying on a surface simultaneously as:
- A 3D BSpline curve (the actual 3D shape).
- A 2D BSpline curve in the surface's parameter space.

The evaluator calls both `D0` in 3D and `Value` in 2D (UV parameter) and feeds both into the 1D / 2D / 3D sub-spaces of `AdvApprox_ApproxAFunction`.

### `Approx_SweepApproximation`

Approximates a sweep surface defined by a **section law** (a family of curves `C(u)` parameterized by the spine parameter `v`):

1. The user provides `Approx_SweepFunction`, implementing `D0(v)` → BSpline section curve at `v`, `D1(v)` → derivative of section.
2. `Approx_SweepApproximation::Perform()` creates an internal `EvaluatorFunction` adaptor.
3. Calls `AdvApprox_ApproxAFunction` to fit the surface in `v` while keeping the `u` polynomial structure from the section curves.
4. The result is a tensor-product BSpline surface `S(u,v)`.

The evaluator can return D0, D1, D2 derivatives in `v`, allowing C0/C1/C2 continuity control.

### `Approx_MCurvesToBSpCurve`

Joins a sequence of `AppParCurves_MultiCurve` (Bezier segments) into a single `AppParCurves_MultiBSpCurve` by:
1. Computing the knot vector from the Bezier segment boundaries.
2. Converting Bezier poles to BSpline poles using degree elevation / pole insertion.

### `Approx_SameParameter`

Given two curves that geometrically coincide but have different parametrizations, computes a reparametrization function (a 1D BSpline map) to make them share the same arc-length parametrization. Uses iterative bisection + Newton refinement.

### `Approx_CurvilinearParameter` and `Approx_CurvlinFunc`

`Approx_CurvlinFunc` wraps a curve or curve-on-surface and provides evaluation of its arc-length function `s(u)`.
`Approx_CurvilinearParameter` fits this arc-length function by a BSpline so the result can be used as a reparametrization map.

### Generic Templates

These `.gxx` files are the generic compute loops instantiated by other packages:

| Template | Instantiated by | Purpose |
|---|---|---|
| `Approx_ComputeLine.gxx` | `BRepApprox`, etc. | Bezier multi-curve fit loop |
| `Approx_ComputeCLine.gxx` | `BRepApprox`, etc. | Continuous-function fit |
| `Approx_BSplComputeLine.gxx` | `BRepApprox`, etc. | BSpline multi-curve fit loop |
| `Approx_FitAndDivide.gxx` | 3D packages | Fit + subdivide for 3D |
| `Approx_FitAndDivide2d.gxx` | 2D packages | Fit + subdivide for 2D |

## Files

### `Approx_Curve3d.hxx / .cxx`
Simplest entry point. Constructor: `(Adaptor3d_Curve, Tol3d, Order, MaxSegments, MaxDegree)`. Result: `Curve()` returns `Geom_BSplineCurve`.

### `Approx_Curve2d.hxx`
2D analogue of `Approx_Curve3d`. Works with `Adaptor2d_Curve2d`.

### `Approx_CurveOnSurface.hxx / .cxx`
Simultaneous 3D + 2D approximation. Constructor takes `Adaptor3d_CurveOnSurface`. Returns both 3D and 2D (PCurve) BSpline curves.

### `Approx_SweepApproximation.hxx / .cxx / .lxx`
Section-law surface fitter. Key method: `Perform(First, Last, Tol3d, BoundTol, Tol2d, TolAngular, Continuity, Degmax, Segmax)`. Internal `Eval()` implements the `AdvApprox_EvaluatorFunction` interface.

### `Approx_SweepFunction.hxx / .cxx`
Abstract base class for the section law. Users override to define the section curve at each spine parameter.

### `Approx_MCurvesToBSpCurve.hxx / .cxx`
Post-processing: merge Bezier segments into BSpline. Used by `Compute` algorithms after fitting.

### `Approx_SameParameter.hxx / .cxx`
Reparametrization alignment. Used in surface intersection and curve-on-surface construction.

### `Approx_CurvilinearParameter.hxx / .cxx`
Arc-length reparametrization. Uses `Approx_CurvlinFunc` as the function to be approximated.

### `Approx_CurvlinFunc.hxx / .cxx`
Arc-length integrator using Gauss quadrature. Provides `D0()` (arc-length at parameter) and `D1()` (speed).

### `Approx_ParametrizationType.hxx`
Enum: `ChordLength`, `Centripetal`, `IsoParametric`.

### `Approx_Status.hxx`
Enum: `Approx_PartialAtBest`, `Approx_NoPointsAdded`, `Approx_NoApproximation`.

### `Approx_FitAndDivide.hxx / _0.cxx` and `FitAndDivide2d.hxx / _0.cxx`
Instantiations of the fit-and-divide generic loop. Used internally by compute classes.

## Summary Table

| File | Role | Key Concept |
|---|---|---|
| `Approx_Curve3d` | 3D curve approximation | Adaptor → BSpline via AdvApprox |
| `Approx_Curve2d` | 2D curve approximation | Adaptor2d → BSpline2d |
| `Approx_CurveOnSurface` | 3D + PCurve simultaneous | Multi-space AdvApprox |
| `Approx_SweepApproximation` | Sweep surface | Section-law → BSpline surface |
| `Approx_SweepFunction` | Section law interface | User-defined section at spine param |
| `Approx_MCurvesToBSpCurve` | Bezier → BSpline join | Segment assembly |
| `Approx_SameParameter` | Reparametrization | Arc-length matching |
| `Approx_CurvilinearParameter` | Arc-length param | BSpline arc-length map |
| `Approx_CurvlinFunc` | Arc-length integrator | Gauss quadrature speed |
| `Approx_ParametrizationType.hxx` | Enum | Chord / Centripetal / Uniform |
| `Approx_Status.hxx` | Enum | Success / partial / failure |
| `Approx_FitAndDivide[2d]` | Generic fit loop | Instantiation of gxx templates |
| `Approx_Compute*.gxx` | Generic fitters | Template algorithms |
