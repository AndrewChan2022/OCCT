# TKTopAlgo / BRepApprox — BRep Topology Fitting

## Architecture

```
BRepApprox_Approx             (main driver: fits intersection lines on BRep faces)
    |
    +-- BRepApprox_ApproxLine  (input: a discrete intersection line with 3D + 2D points)
    |
    +-- Generated fitting classes (same pattern as AppDef, instantiated for BRepApprox types):
         BRepApprox_MyBSplGradientOfTheComputeLineOfApprox
         BRepApprox_MyGradientOfTheComputeLineBezierOfApprox
         BRepApprox_MyGradientbisOfTheComputeLineOfApprox
         BRepApprox_BSpParLeastSquare*
         BRepApprox_BSpParFunction*
         BRepApprox_BSpGradient_BFGS*
         BRepApprox_Gradient_BFGS*
         (... and 2D variants)
```

## Purpose

`BRepApprox` (inside `TKTopAlgo`) handles the **approximation of curve intersections** on BRep solid geometry.

When two BRep faces intersect, the intersection result is a set of 3D curve points plus their images in the 2D parameter spaces of each face (PCurves). `BRepApprox` fits all three representations simultaneously:
- 1 × 3D BSpline curve.
- 2 × 2D BSpline PCurves (one per face).

This multi-space simultaneous fit ensures geometric consistency: the 3D curve and the two PCurves are all fitted together so that their errors are balanced.

## Principle / Theory / Algorithm

The algorithm is identical to `AppDef_BSplineCompute` / `AppParCurves` but instantiated with BRep-specific types:

1. **Input**: `BRepApprox_ApproxLine` stores the discrete intersection points as a MultiLine-like structure containing 3D points and 2D parameter points on each face.

2. **Fitting**: `BRepApprox_Approx` runs the same BFGS-parameterization + BSpline least-squares loop as `AppDef_BSplineCompute`, fitting:
   - 1 sub-space of dimension 3 (the 3D curve).
   - 2 sub-spaces of dimension 2 (the two PCurves).

3. **Tolerance control**: tolerances are specified independently for the 3D curve (`Tol3d`) and the 2D PCurves (`Tol2d`). The 2D tolerance is typically derived from the surface resolution at the intersection.

4. **Subdivision**: if the fitting error exceeds tolerance, the intersection line is subdivided and the two sub-segments are fitted independently.

5. **Output**: multiple `AppParCurves_MultiBSpCurve` objects, one per segment. The caller assembles them into a full intersection result.

### Why separate from `AppDef`?

`AppDef` provides `AppDef_MultiLine` / `AppDef_MyLineTool` as concrete types. `BRepApprox` provides its own `BRepApprox_TheMultiLine` (hidden behind generated code) that reads intersection data differently. The `.gxx` template mechanism allows full code reuse with different concrete line types without virtual dispatch overhead.

## Files

### `BRepApprox_ApproxLine.hxx / .cxx`
**Input data container.**
Stores:
- `NbPoints()`: number of intersection points.
- `Point(i)`: the 3D point at sample `i`.
- `Pnt2dOnS1(i)`, `Pnt2dOnS2(i)`: 2D parameter-space images on each face.
- `HasTangency()`, `Tangent(i)`, `TangentOnS1(i)`, `TangentOnS2(i)`: optional tangent vectors.

### `BRepApprox_Approx.hxx / _0.cxx`
**Main driver.**
Instantiates `Approx_BSplComputeLine.gxx` (via the `_0.cxx` include) with `BRepApprox` types. Key method: `Perform(Line, Tol3d, Tol2d)`. Returns:
- `Value()`: `AppParCurves_MultiBSpCurve` with 1 × 3D + 2 × 2D curves.
- `Error(tol3d, tol2d)`: achieved tolerances.

### Generated `BRepApprox_*_0.cxx` files
One file per `.gxx` generic template instantiation. No algorithmic content — just template expansion hooks. This mirrors the pattern of `AppDef_*_0.cxx`.

### Gradient / BFGS classes (`BRepApprox_*Gradient*`, `BRepApprox_*BFGS*`)
Instantiations of `AppParCurves_Gradient.gxx`, `AppParCurves_BSpGradient.gxx`, and the corresponding BFGS wrappers, specialised for BRepApprox line types.

### Least-squares classes (`BRepApprox_*LeastSquare*`)
Instantiations of `AppParCurves_LeastSquare.gxx` for Bezier and BSpline, specialised for BRepApprox line types.

## Summary Table

| File | Role | Key Concept |
|---|---|---|
| `BRepApprox_ApproxLine` | Intersection line input | 3D + 2D PCurve sample data |
| `BRepApprox_Approx` | Main driver | Multi-space BSpline fit |
| `BRepApprox_*LeastSquare*` | LS solver instances | Normal equations for BRep data |
| `BRepApprox_*Function*` | Objective instances | E(params) for BRep types |
| `BRepApprox_*Gradient*` | Gradient instances | ∂E/∂u for BRep types |
| `BRepApprox_*BFGS*` | BFGS optimizer | Parameter update for BRep |
| (all `_0.cxx` files) | Template instantiation | No code, just `#include` |
