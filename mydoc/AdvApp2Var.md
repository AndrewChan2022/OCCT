# AdvApp2Var — Adaptive 2D Surface Approximation

## Architecture

```
AdvApp2Var_ApproxAFunc2Var   (top-level driver)
    |
    +-- AdvApp2Var_Context       (approximation parameters / precision codes)
    |
    +-- AdvApp2Var_Network       (patch grid: collection of Patches)
    |     +-- AdvApp2Var_Patch   (one rectangular domain [Ui,Ui+1] x [Vj,Vj+1])
    |
    +-- AdvApp2Var_Framework     (boundary iso-curves and node constraints)
    |     +-- AdvApp2Var_Iso     (one parametric iso-line approximation)
    |     +-- AdvApp2Var_Node    (corner-point constraint values)
    |
    +-- AdvApp2Var_Criterion     (optional user convergence criterion)
    |
    +-- AdvApp2Var_MathBase      (low-level Fortran-translated math routines)
    +-- AdvApp2Var_SysBase       (Fortran I/O stubs)
    +-- AdvApp2Var_Data_*.pxx    (pre-computed Jacobi/quadrature tables)
```

`AdvApprox_Cutting` (from `AdvApprox`) is reused for knot insertion decisions in U and V.

## Purpose

`AdvApp2Var` approximates a function `F(U,V) : ℝ² → ℝ³` (and optionally 1D/2D sub-spaces) by a **tensor-product BSpline surface**, within per-dimension tolerances and boundary-continuity tolerances.

It extends `AdvApprox` (1D) to 2 parameters by:
1. Fixing V isoparametric lines and approximating them as 1D curves in U.
2. Assembling the results patch-by-patch.
3. Adaptively subdividing in U and V.

## Principle / Theory / Algorithm

### 1. Tensor-Product Decomposition

The 2D problem is reduced to repeated 1D problems along iso-lines:

```
F(U, V) is first discretised along a grid of V-values (V-iso lines).
For each V-iso line, the 1D AdvApprox machinery fits F(·, V) in U.
The U-polynomial coefficients are then fit in V.
```

### 2. Patch Grid

The domain `[U₀,U₁] × [V₀,V₁]` is divided into a `Network` of rectangular `Patch` objects. Each patch stores polynomial coefficients in Jacobi basis for its sub-domain.

### 3. Iso-Curve Fitting (`AdvApp2Var_Iso`)

An iso-curve (fixed U or fixed V) is approximated using the same Gauss-Jacobi quadrature as `AdvApprox_SimpleApprox`. Results (Jacobi coefficients) are stored on the `Framework` boundary structure.

### 4. Node Constraints (`AdvApp2Var_Node`)

At patch corners, function values and derivatives computed from the iso-curve fits are stored as point constraints to ensure C0/C1/C2 continuity across patch boundaries.

### 5. Adaptive Subdivision

Convergence is tested per-patch:
- If a patch does not meet the tolerance (or fails the optional `AdvApp2Var_Criterion`), it is subdivided using `AdvApprox_Cutting` in U or V.
- The subdivision preference (`GeomAbs_IsoType FavorIso`) biases the algorithm toward inserting knots in U or V first.

### 6. Precision Code

`PrecisCode ∈ {1, 2, 3}` selects the number of Gauss points used:
- 1 → fast, average precision
- 2 → balanced
- 3 → slow, high precision

### 7. BSpline Conversion

After all patches converge, `ConvertBS()` assembles the patch polynomial arrays into a global `Geom_BSplineSurface` by concatenating knot vectors and converting Jacobi coefficients to B-Spline poles.

## Files

### `AdvApp2Var_EvaluatorFunc2Var.hxx`
Abstract evaluator for `F(U,V)`. The user sub-classes this and implements `operator()` with signature matching the expected parameter ordering.

### `AdvApp2Var_Context.hxx / .cxx`
Stores approximation parameters: degree limits, continuity orders, precision code, parameter bounds. Acts as a configuration object passed to all sub-algorithms.

### `AdvApp2Var_Criterion.hxx / .cxx`
Optional user-defined convergence criterion. Override `IsSatisfied(Patch)` to add problem-specific error metrics beyond the default L-infinity bound.

### `AdvApp2Var_CriterionRepartition.hxx` / `CriterionType.hxx`
Enumerations controlling how the criterion guides subdivision (which direction, how to repartition).

### `AdvApp2Var_Network.hxx / .cxx`
Container for the 2D array of `Patch` objects. Manages insertion and iteration over patches during adaptive refinement.

### `AdvApp2Var_Patch.hxx / .cxx`
Represents one `[Ui,Ui+1] × [Vj,Vj+1]` domain.
- Stores polynomial coefficients (Jacobi basis) once the patch is computed.
- Tracks whether it is discretised and whether it has converged.
- `ComputeApprox(Context, Framework, Func)`: performs the 1D fits on iso-lines within this patch.

### `AdvApp2Var_Framework.hxx / .cxx`
Holds the collection of `Iso` and `Node` objects forming the boundaries and corners of all patches. Provides access by index.

### `AdvApp2Var_Iso.hxx / .cxx`
Approximates one parametric iso-line. Calls Fortran-translated `MathBase` routines to perform Gauss-Jacobi quadrature on the iso-line samples.

### `AdvApp2Var_Node.hxx / .cxx`
Stores function value and derivatives at a patch corner point. Used to enforce continuity when assembling the global surface.

### `AdvApp2Var_MathBase.hxx / .cxx`
Low-level mathematical routines translated from Fortran (originally from the MIRIAD library). Includes Gauss point computation, matrix operations for polynomial conversion, and Jacobi basis evaluation.

### `AdvApp2Var_SysBase.hxx / .cxx`
Portability stubs for Fortran-style I/O and error handling used by `MathBase`.

### `AdvApp2Var_Data_*.pxx`
Pre-computed numerical tables (Gauss weights, Jacobi polynomial values, knot sequences) embedded as C++ data files. Generated once and compiled in.

### `AdvApp2Var_ApproxF2var.hxx / .cxx`
Internal helper calling the Fortran-translated core approximation routines for a single patch. Bridge between C++ classes and the Fortran-origin numerical core.

### `AdvApp2Var_ApproxAFunc2Var.hxx / .cxx / .lxx`
**Top-level driver.**
- Constructor wires together `Context`, `Network`, `Framework`, evaluator, and cutting tools.
- `Perform()` / `ComputePatches()` / `ComputeConstraints()`: orchestrate the adaptive loop.
- `Compute3DErrors()`: post-processing error measurement.
- `ConvertBS()`: assemble the final `Geom_BSplineSurface`.
- Result: `Surface(Index)` returns the fitted `Geom_BSplineSurface`.

## Summary Table

| File | Role | Key Concept |
|---|---|---|
| `AdvApp2Var_EvaluatorFunc2Var.hxx` | Abstract 2D evaluator | User callback F(U,V) |
| `AdvApp2Var_Context` | Approximation parameters | Degree, continuity, precision |
| `AdvApp2Var_Criterion` | Optional convergence test | User-defined quality metric |
| `AdvApp2Var_Network` | Patch grid container | 2D array of patches |
| `AdvApp2Var_Patch` | One rectangular domain | Per-patch Jacobi coefficients |
| `AdvApp2Var_Framework` | Boundary iso + node store | Continuity data structure |
| `AdvApp2Var_Iso` | Iso-line fit | 1D Gauss-Jacobi on iso |
| `AdvApp2Var_Node` | Corner constraint | F value + derivatives at corner |
| `AdvApp2Var_ApproxF2var` | Internal core bridge | Fortran-C++ interface |
| `AdvApp2Var_MathBase` | Low-level math | Gauss points, polynomial tables |
| `AdvApp2Var_SysBase` | Fortran portability stubs | I/O and error stubs |
| `AdvApp2Var_Data_*.pxx` | Pre-computed tables | Numerical constants |
| `AdvApp2Var_ApproxAFunc2Var` | Top-level driver | Adaptive 2D approximation loop |
