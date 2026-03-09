# AdvApprox — Adaptive 1D Function Approximation

## Architecture

```
AdvApprox_ApproxAFunction   (top-level driver)
    |
    +-- AdvApprox_SimpleApprox   (single-segment fit: Hermite + Jacobi)
    |
    +-- AdvApprox_Cutting        (abstract: how to bisect an interval)
         +-- AdvApprox_DichoCutting   (binary midpoint)
         +-- AdvApprox_PrefAndRec     (preferred + recommended cut points)
         +-- AdvApprox_PrefCutting    (preferred cut points)
```

`AdvApprox_EvaluatorFunction` is the abstract callback that the caller must supply to evaluate the function being approximated.

## Purpose

`AdvApprox` approximates one or more scalar/vector functions of one parameter `t ∈ [First, Last]` by a **piecewise polynomial BSpline** within user-specified tolerances.
It supports simultaneous fitting of 1D scalar, 2D vector, and 3D vector sub-spaces.

The result is a single knot vector shared by all sub-spaces, with one common polynomial degree.
This module is the backbone used by `AdvApp2Var` (surface approximation) and `Approx_SweepApproximation`.

---

## Single-Segment Algorithm: Hermite + Jacobi Decomposition

This is the core mathematical idea, implemented in `AdvApprox_SimpleApprox::Perform()`.

Each segment fits the function on a normalized interval `[-1, 1]` using a **two-part decomposition**:

```
P(t) = R(t)  +  W(t) · Q(t)
        ↑              ↑
   Hermite part    Jacobi part
   (endpoints)     (interior)
```

### Part 1 — R(t): Hermite Interpolation at Endpoints

The function and its derivatives up to order `NivConstr` are evaluated **exactly** at both ends of the segment:
- At `t = First`: `f(First)`, `f'(First)`, ..., `f^(NivConstr)(First)`
- At `t = Last`:  `f(Last)`,  `f'(Last)`,  ..., `f^(NivConstr)(Last)`

This gives `2*(NivConstr+1)` values in total (matching conditions at both ends).

`PLib::HermiteInterpolate` constructs the unique polynomial `R(t)` of degree `DegreeR = 2*NivConstr + 1` that satisfies all these endpoint conditions **exactly**.

```
NivConstr = 0  →  DegreeR = 1  (linear, C0 endpoints)
NivConstr = 1  →  DegreeR = 3  (cubic,  C1 endpoints)
NivConstr = 2  →  DegreeR = 5  (quintic, C2 endpoints)
```

The Hermite polynomial `R(t)` is stored in `myFirstConstr` / `myLastConstr` (the endpoint data) and `myCoeff` (the resulting polynomial coefficients).

**Why this is essential for inter-segment continuity:** because the endpoint constraints come directly from evaluating the *original function*, adjacent segments that share a boundary point will both produce Hermite terms that exactly reproduce the function value (and derivatives) at that point. This makes the joint automatically match — not approximately, but exactly up to `NivConstr` order — regardless of the interior fit quality.

### Part 2 — W(t) · Q(t): Jacobi Polynomial for the Interior Residual

The residual `F(t) - R(t)` is zero at the endpoints (by construction of R). It is approximated by a polynomial `Q(t)` multiplied by the **Jacobi weight function**:

```
W(t) = (1 - t²)^(NivConstr + 1)
```

`W(t)` is exactly zero at `t = ±1` with multiplicity `NivConstr+1`, so `W(t)·Q(t)` cannot disturb any derivative up to order `NivConstr` at either endpoint. Adding it to `R(t)` preserves the Hermite constraints perfectly.

`Q(t)` is expanded in the **Jacobi orthogonal polynomial basis** `{J₀, J₁, …, J_DegreeQ}`, which is orthogonal with weight `W(t)` on `[-1, 1]`:

```
Q(t) = Σₖ cₖ Jₖ(t)
```

The coefficients `cₖ` are computed by **Gauss–Jacobi quadrature** — evaluating the function at Gauss points `{tᵢ}` and taking weighted inner products:

```
cₖ = Σᵢ wᵢₖ · [F(tᵢ) - R(tᵢ)]
```

where `wᵢₖ` are the pre-tabulated Gauss–Jacobi weights from `PLib_JacobiPolynomial`. This replaces a linear system solve with a simple dot product.

The total polynomial degree is:
```
WorkDegree = DegreeR + DegreeQ + 2*(NivConstr+1)
           = (2*NivConstr+1) + DegreeQ + (2*NivConstr+2)
```

### Degree Reduction and Error Estimation

After computing Jacobi coefficients `{cₖ}` at full `WorkDegree`, `PLib_JacobiPolynomial::ReduceDegree()` finds the **minimum degree** `NewDegree` such that the truncated series still satisfies the tolerance. The error introduced by dropping high-degree Jacobi terms is measured via `MaxError()` and `AverageError()`.

A minimum degree floor is enforced when subdivision has occurred:
```cpp
if (isCut && TheDeg < 2 * ContinuityOrder + 1)
    TheDeg = 2 * ContinuityOrder + 1;
```
This prevents the polynomial from being too low-degree to faithfully encode the Hermite constraints at both ends.

### Summary of `SimpleApprox::Perform()` Steps

```
1. Evaluate F and its derivatives at First and Last  → myFirstConstr, myLastConstr
2. PLib::HermiteInterpolate(...)                      → R(t) in myCoeff[0..DegreeR]
3. Evaluate F at Gauss points {tᵢ}, compute F(tᵢ) - R(tᵢ)
   (using mySomTab / myDifTab for symmetric + antisymmetric parts)
4. Gauss-Jacobi quadrature                            → Jacobi coefficients cₖ
   stored in myCoeff[DegreeR+1 .. WorkDegree]
5. ReduceDegree()                                     → find minimal NewDegree
6. MaxError / AverageError                            → error bounds
```

---

## Multi-Segment Algorithm: Adaptive Subdivision Loop

`AdvApprox_ApproxAFunction::Approximation()` drives the full approximation:

### Interval Stack

The driver maintains a stack of sub-intervals (stored as `TABINT[]`). Initially `TABINT = {First, Last}`.

### Main Loop

```
while stack is not exhausted:
    1. Pop current interval [a, b]
    2. Call SimpleApprox.Perform([a, b])  → one polynomial segment
    3. If error ≤ tolerance for all sub-spaces:
           accept segment, store coefficients
    4. Else if nSegments < MaxSeg:
           cut [a, b] into [a, c] and [c, b] via CutTool
           push both back onto the stack
    5. Else (max segments reached):
           accept anyway (HasResult = true, IsDone = false)
```

The "stack" is implemented as a flat sorted array `TABINT[0..NUPIL]` with the current interval always at `TABINT[NumCurves]`.

### After the Loop: `PrepareConvert()` — Continuity Check at Joints

After all segments are fitted, `PrepareConvert()` examines each **inter-segment joint** between segment `i` and segment `i+1`. At joint `k = TABINT[i]`:

1. **Evaluate** segment `i` at its right endpoint and segment `i+1` at its left endpoint, including all derivatives up to `ContinuityOrder`. This uses `PLib::EvalPolynomial` on the raw polynomial coefficient arrays.

2. **Scale the derivatives** by the chain-rule factor to convert from normalized `[-1,1]` to physical parameter space:
   ```
   facteurᵢ = ((bᵢ - aᵢ) / 2)^order
   ```
   where `[aᵢ, bᵢ]` is the physical interval of segment `i`.

3. **Compare** the scaled derivatives from both sides. If the difference is within a small fraction of the local tolerance, the joint achieves continuity order `iordre`. The residual error is accumulated into `ErrorMax`.

4. **Record** the achieved continuity in `TabContinuity[joint]` — this is the actual delivered continuity at that knot, which may be anywhere from C0 to C(ContinuityOrder).

### Why the Joint is (Usually) Exact

Because each `SimpleApprox` evaluated `F` and its derivatives **exactly** at its endpoints for the Hermite part `R(t)`, both segments independently produce a polynomial whose endpoint values match the true function. If the function is at least `C_NivConstr` continuous at the joint, the two-sided Hermite values agree and the joint is automatically `C_NivConstr` up to floating-point precision.

If `F` has a genuine discontinuity (e.g., a kink), the Hermite values on the two sides will differ, `PrepareConvert` will detect this, and the joint will be recorded as C(-1) (a break), leading to a knot with full multiplicity.

### `Convert_CompPolynomialToPoles()` — Building the Final BSpline

`Convert_CompPolynomialToPoles` takes:
- All polynomial segments (Taylor coefficients on normalized `[-1,1]`).
- The physical knot positions `IntervalsArray`.
- `TabContinuity[i]` — the continuity achieved at each joint.

It produces a single BSpline by:
1. **Setting knot multiplicities** from the continuity: at joint `i`, multiplicity = `degree - TabContinuity[i]`.
   - C2 joint → multiplicity = `degree - 2`
   - C1 joint → multiplicity = `degree - 1`
   - C0 joint → multiplicity = `degree`   (a full break, curve passes through)
2. **Converting polynomial coefficients** to BSpline control poles using `BSplCLib` routines (knot insertion / Greville / direct formula).

The final output `myKnots`, `myMults`, `myDegree`, and `my[1/2/3]DPoles` represents the complete piecewise BSpline.

### Complete Data Flow

```
User function F(t)
      |
      | (evaluated at endpoints)       (evaluated at Gauss points)
      ↓                                        ↓
 Hermite R(t)  ──────────── F(tᵢ) - R(tᵢ)  ──→  Jacobi Q(t)
      |                                               |
      └──────── P(t) = R(t) + W(t)·Q(t) ─────────────┘
                          |
                   [for each segment]
                          |
              PrepareConvert: check joint derivatives
                          |
                   TabContinuity[i]
                          |
         Convert_CompPolynomialToPoles
                          |
              BSpline (knots + poles)
```

---

## Files

### `AdvApprox_EvaluatorFunction.hxx`
**Abstract callback interface.**
The user sub-classes this and implements `operator()(Dimension, FirstLast, Parameter, DerivativeOrder, Result, ErrorCode)`.
- `DerivativeOrder = 0` → position values.
- `DerivativeOrder = 1` → first derivatives (scaled by `(Last-First)/2`).
- `DerivativeOrder = 2` → second derivatives (scaled by `((Last-First)/2)²`).

Result array layout: `[1D subspaces | 2D subspaces (×2) | 3D subspaces (×3)]`.

### `AdvApprox_SimpleApprox.hxx / .cxx`
**Single-segment Hermite + Jacobi fit.**

Key members:
- `myNivConstr`: continuity order (0/1/2), governs Hermite degree and Jacobi weight.
- `myFirstConstr` / `myLastConstr`: `[dim × (NivConstr+1)]` arrays of endpoint values and derivatives.
- `myCoeff`: full polynomial coefficients `[0..WorkDegree × Dim]` combining Hermite + Jacobi parts.
- `myTabPoints`: Gauss abscissae from `PLib_JacobiPolynomial::Points()`.
- `myTabWeights`: Gauss weights from `PLib_JacobiPolynomial::Weights()`.
- `mySomTab`, `myDifTab`: symmetric and antisymmetric halves of `F(tᵢ) - R(tᵢ)` used to separate even/odd Jacobi coefficients.

Key flow in `Perform()`:
1. Evaluate `F` at `First` and `Last` for orders `0..NivConstr` → store in `myFirstConstr`, `myLastConstr`.
2. `PLib::HermiteInterpolate` on normalized `[-1,1]` → Hermite polynomial `R(t)` in `myCoeff`.
3. Evaluate `F` at Gauss points; subtract `R(tᵢ)` using `PLib::EvalPolynomial`.
4. Weighted sum over Gauss points → Jacobi coefficients appended to `myCoeff`.
5. `ReduceDegree()` → find minimum degree satisfying tolerance.
6. `MaxError()` / `AverageError()` → error bound per sub-space.

### `AdvApprox_ApproxAFunction.hxx / .cxx / .lxx`
**Multi-segment adaptive driver.**

Key private methods:
- `Approximation()`: the main subdivision loop (interval stack, SimpleApprox calls, cutting).
- `PrepareConvert()`: joint-by-joint continuity checker; sets `TabContinuity`.
- Uses `Convert_CompPolynomialToPoles` to assemble the final BSpline.

Two constructors: default `DichoCutting`, or user-supplied `CutTool`.
Static `Approximation(...)` allows calling the loop without constructing an object.

### `AdvApprox_Cutting.hxx`
**Abstract cut strategy.** Single pure-virtual `Value(a, b, cuttingvalue)`.

### `AdvApprox_DichoCutting.hxx / .cxx`
**Binary midpoint.** Always returns `(a + b) / 2`.

### `AdvApprox_PrefCutting.hxx / .cxx`
**Preferred-parameter list.** Finds the listed parameter nearest the midpoint of `[a, b]`.

### `AdvApprox_PrefAndRec.hxx / .cxx`
**Preferred + recommended.** Preferred list has priority; recommended list is the fallback. Useful when discontinuities are known in advance — put them in the preferred list so the knot lands exactly there.

---

## Summary Table

| File | Role | Key Concept |
|---|---|---|
| `AdvApprox_EvaluatorFunction.hxx` | Abstract evaluator callback | User-implemented D0/D1/D2 oracle |
| `AdvApprox_SimpleApprox` | Single-segment fit | **Hermite** (endpoints) + **Jacobi** (interior) |
| `AdvApprox_ApproxAFunction` | Adaptive multi-segment driver | Interval stack, subdivision, `PrepareConvert`, `Convert_CompPolynomialToPoles` |
| `AdvApprox_Cutting` | Abstract cut strategy | Interval bisection policy |
| `AdvApprox_DichoCutting` | Midpoint cut | `c = (a+b)/2` |
| `AdvApprox_PrefCutting` | Preferred-point cut | Nearest preferred parameter |
| `AdvApprox_PrefAndRec` | Preferred + recommended cut | Priority-list, useful at known discontinuities |

---

## Key Design Insight

The split `P = R + W·Q` is elegant because the two parts have completely separate responsibilities:

- **R(t)** (Hermite): handles **inter-segment continuity**. It is the only part that touches endpoint derivatives and it is computed exactly, not approximately.
- **W(t)·Q(t)** (Jacobi): handles **interior approximation quality**. It is zero at the endpoints so it cannot break continuity, and it uses an orthogonal basis so Gauss quadrature gives the optimal coefficients without a matrix solve.

This separation means continuity at joints is achieved structurally (by construction), not by tuning tolerances — and the interior quality is handled independently by the Jacobi expansion.
