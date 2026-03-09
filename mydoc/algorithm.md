# Key Algorithms — OCCT Fitting & Approximation

## 1. Continuous Rebuild (AdvApprox path)

Input: an evaluable function `F(t)` (or `F(u,v)` for surfaces).

### 1a. Single Segment Fitting — Hermite + Jacobi Decomposition

**Core idea**: split the polynomial into two parts with separated responsibilities.

```
P(t) = R(t)  +  W(t) * Q(t)
        |              |
   Hermite part    Jacobi part
   (endpoints)     (interior)
```

**Step 1 — Hermite interpolation R(t)**

Evaluate `F` and its derivatives at both endpoints up to order `NivConstr`:

```
NivConstr = 0  =>  R is degree 1  (linear,  matches C0 at endpoints)
NivConstr = 1  =>  R is degree 3  (cubic,   matches C1 at endpoints)
NivConstr = 2  =>  R is degree 5  (quintic, matches C2 at endpoints)
```

R(t) reproduces the function exactly at endpoints. This is what guarantees
inter-segment continuity — adjacent segments both evaluate the same function
at the shared boundary, so they agree automatically up to `NivConstr` order.

**Step 2 — Jacobi correction W(t)*Q(t)**

The residual `F(t) - R(t)` is zero at the endpoints (by construction).
Multiply by the Jacobi weight:

```
W(t) = (1 - t^2)^(NivConstr + 1)
```

W(t) has `NivConstr+1` zeros at t = +/-1, so the product `W*Q` cannot
disturb any derivative at the endpoints up to order NivConstr. The Hermite
continuity is preserved perfectly.

Q(t) is expanded in the **Jacobi orthogonal polynomial basis** `{Jk}`:

```
Q(t) = sum_k  c_k * J_k(t)
```

Coefficients `c_k` are computed by **Gauss-Jacobi quadrature** — just
weighted dot products at Gauss points. No matrix solve needed.

**Step 3 — Degree reduction**

After computing all Jacobi coefficients up to `WorkDegree`, drop
high-degree terms until the truncation error exceeds tolerance.
This gives the minimum polynomial degree per segment.

**Key file**: `AdvApprox_SimpleApprox` — implements the full Hermite+Jacobi
single-segment pipeline.

### 1b. Adaptive Segment Subdivision

The driver `AdvApprox_ApproxAFunction` manages a stack of intervals:

```
TABINT = [First, Last]           <- initial single interval

while stack not empty:
    pop [a, b]
    fit with SimpleApprox
    if error <= tolerance:
        accept segment
    else if nSegments < MaxSeg:
        cut [a, b] -> [a, c] + [c, b]    <- via Cutting strategy
        push both halves
    else:
        accept with degraded tolerance
```

**Cutting strategies** (where to place the cut point `c`):

| Strategy             | Cut rule                                      |
|----------------------|-----------------------------------------------|
| `DichoCutting`       | Midpoint: `c = (a+b)/2`                       |
| `PrefCutting`        | Nearest preferred parameter to midpoint        |
| `PrefAndRec`         | Preferred list first, recommended as fallback  |

The preferred-parameter strategies are useful when known discontinuities
(kinks, C0/C1 breaks) should land exactly on knot positions.

**After fitting all segments — joint verification (`PrepareConvert`)**:

At each inter-segment joint:
1. Evaluate both segments' derivatives at the shared boundary.
2. Scale by chain-rule factor `((b-a)/2)^order`.
3. Compare. Record achieved continuity: C0, C1, C2, or break.

The joint is usually exact because both segments' Hermite parts evaluated
the same function at the same point. Discrepancy only arises if the function
has a genuine discontinuity.

**Final assembly (`Convert_CompPolynomialToPoles`)**:

- Set knot multiplicity from continuity: `mult = degree - continuity_order`.
- Convert polynomial coefficients to BSpline control poles.

**Key files**:
- `AdvApprox_ApproxAFunction` — adaptive loop driver
- `AdvApprox_Cutting` / `DichoCutting` / `PrefCutting` / `PrefAndRec` — cut strategies

### 1c. Surface Extension (AdvApp2Var)

Extends the 1D algorithm to 2D by tensor-product decomposition:

```
1. Divide [U0,U1] x [V0,V1] into rectangular patches
2. For each patch, fix V iso-lines and fit 1D in U (same Hermite+Jacobi)
3. Assemble patch Jacobi coefficients
4. Check convergence per patch; if fail, subdivide in U or V
5. Convert all patches to a global BSpline surface
```

Corner `Node` constraints and boundary `Iso` fits enforce continuity
across patch boundaries, mirroring the 1D Hermite endpoint strategy.

**Key file**: `AdvApp2Var_ApproxAFunc2Var` — adaptive 2D driver

---

## 2. Discrete Point Fitting (AppDef path)

Input: discrete data points `{P_1, ..., P_n}` with optional tangent/curvature.

### 2a. BFGS Parameter Optimization

The parameters `{u_i}` assigned to each point are unknown. The algorithm
alternates between two steps:

```
1. Initialize {u_i} by chord-length / centripetal / uniform

2. Outer loop (BFGS):
   a. Fix {u_i}, solve least-squares for BSpline poles Q:
      (A^T A) Q = A^T P        <- normal equations
      where A is the BSpline collocation matrix at {u_i}

   b. Evaluate objective:
      E = sum_i || C(u_i) - P_i ||^2

   c. Compute gradient analytically:
      dE/du_i = 2 * (C(u_i) - P_i) . C'(u_i)

   d. BFGS update: adjust {u_i} using approximate inverse Hessian

3. Repeat until E < tol^2 or max iterations
```

BFGS only changes the parameters. Poles are re-solved by linear
least-squares at each iteration. This separates the nonlinear part
(parametrization) from the linear part (pole computation).

**Key files**:
- `AppParCurves_LeastSquare.gxx` — normal equations solver
- `AppParCurves_Function.gxx` / `BSpFunction.gxx` — objective E(u)
- `AppParCurves_Gradient.gxx` / `BSpGradient.gxx` — analytic gradient
- `AppDef_BSplineCompute` — BSpline fitter driver
- `AppDef_Compute` — Bezier fitter driver

### 2b. Adaptive Segment Subdivision (AppDef's cutting)

AppDef's cutting is simpler than AdvApprox — no Hermite, no Jacobi:

```
if error > tolerance AND cutting enabled:
    find the point with maximum error
    split the point set into two subsets at that point
    run the full BFGS fit on each subset independently
```

Each subset gets its own BSpline segment. The results are joined via
`Approx_MCurvesToBSpCurve` (Bezier segments -> single BSpline).

### 2c. Variational / Fair Curve (AppDef_Variational)

Adds regularization energy terms to make the curve smooth:

```
E = w1 * E_data  +  w2 * E_length  +  w3 * E_curvature  +  w4 * E_torsion

where:
  E_data      = sum || C(u_i) - P_i ||^2        (fidelity)
  E_length    = integral || C'(t) ||^2 dt        (1st derivative energy)
  E_curvature = integral || C''(t) ||^2 dt       (2nd derivative energy)
  E_torsion   = integral || C'''(t) ||^2 dt      (3rd derivative energy)
```

For fixed knot vector, this is a linear problem in the BSpline poles.
Solved via FEmTool (finite-element stiffness matrix assembly):

1. Assemble stiffness matrix for each energy term over Hermite elements.
2. Add data-fidelity matrix.
3. Solve the linear system for optimal poles.
4. If fidelity error too large, increase `w1` and re-solve.

**Key file**: `AppDef_Variational`

---

## 3. AppCont — Lightweight Continuous-Function Path

A simpler alternative to AdvApprox for continuous functions:

```
1. Sample F(t) at Gauss-Legendre quadrature points
2. Build Bernstein collocation matrix B
3. Solve normal equations: (B^T B) Q = B^T P
4. If error > tol, subdivide and repeat per segment
```

No Hermite/Jacobi decomposition, no BFGS. Just direct least-squares
with pre-computed `B^T B` matrices.

Used by `Approx_FitAndDivide` / `Approx_FitAndDivide2d`.

**Key file**: `AppCont_LeastSquare`

---

## Algorithm Comparison Table


| Aspect               | AdvApprox (rebuild)            | AppDef (discrete)              | AppCont (lightweight)      |
|----------------------|--------------------------------|--------------------------------|----------------------------|
| **Input**            | Continuous function F(t)       | Discrete points {P_i}          | Continuous function F(t)   |
| **Single segment**   | Hermite + Jacobi decomposition | Least-squares normal equations | Gauss-Legendre + LS        |
| **Parameter source** | Gauss quadrature points        | BFGS-optimized {u_i}           | Gauss-Legendre points      |
| **No-matrix trick**  | Gauss-Jacobi dot product       | —                              | Pre-computed B^T B         |
| **Continuity**       | Hermite endpoints (structural) | Knot multiplicity              | Endpoint constraints       |
| **Subdivision**      | Adaptive bisection + Cutting   | Split at max-error point       | Fit-and-divide             |
| **Surface**          | AdvApp2Var (tensor product)    | —                              | —                          |
| **Smoothing**        | —                              | Variational energy terms       | —                          |
| **Typical caller**   | Approx_Curve3d, SweepApprox    | AppDef_BSplineCompute          | Approx_FitAndDivide        |
