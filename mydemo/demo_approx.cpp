// demo_approx.cpp
// Demonstrates 4 approximation (rebuild/fitting) use-cases using OCCT.
//
// Build (example, adjust paths):
//   g++ demo_approx.cpp -I<OCCT_INC> -L<OCCT_LIB> \
//       -lTKernel -lTKMath -lTKG3d -lTKGeomBase -lTKGeomAlgo -o demo_approx
//
// Demos:
//   5. Single curve rebuild (sample D0 + D1 tangent + D2 curvature)
//   6. Curve chain rebuild (multiple segments, matching at joints)
//   7. Single curve rebuild with prescribed discontinuity parameters
//   8. Surface rebuild (sample D0 + first/second derivatives)

#include <Approx_Curve3d.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_TrimmedCurve.hxx>

#include <AdvApprox_ApproxAFunction.hxx>
#include <AdvApprox_DichoCutting.hxx>
#include <AdvApprox_PrefAndRec.hxx>
#include <AdvApprox_EvaluatorFunction.hxx>

#include <AdvApp2Var_ApproxAFunc2Var.hxx>
#include <AdvApp2Var_EvaluatorFunc2Var.hxx>
#include <AdvApprox_Cutting.hxx>

#include <GeomAbs_Shape.hxx>
#include <GeomAbs_IsoType.hxx>
#include <NCollection_HArray1.hxx>
#include <NCollection_HArray2.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------
static void PrintCurveInfo(const Handle(Geom_BSplineCurve)& C, const char* label)
{
  if (C.IsNull()) { std::cout << label << ": FAILED\n"; return; }
  std::cout << label
            << ": degree=" << C->Degree()
            << "  nbPoles=" << C->NbPoles()
            << "  nbKnots=" << C->NbKnots()
            << "\n";
}

// ===========================================================================
// Demo 5: Single curve rebuild — D0 + D1 + D2 sampling via Approx_Curve3d
// ===========================================================================
// Source curve: a 3D helix arc defined as a parametric function.
// We wrap it as an Adaptor3d_Curve and call Approx_Curve3d.
// The Approx_Curve3d internally uses AdvApprox_ApproxAFunction which samples
// D0, D1, D2 of the adaptor as needed to satisfy the requested continuity.
//
// Here we fit with C2 continuity so D0+D1+D2 are all used to ensure smooth
// patch junctions.
// ===========================================================================

// A simple helix adaptor: P(t) = (R*cos t, R*sin t, H*t/(2π))
class HelixAdaptor : public Adaptor3d_Curve
{
public:
  HelixAdaptor(double R, double H, double t0, double t1)
    : myR(R), myH(H), myT0(t0), myT1(t1) {}

  double FirstParameter() const override { return myT0; }
  double LastParameter()  const override { return myT1; }

  void D0(const double t, gp_Pnt& P) const override {
    P.SetCoord(myR * std::cos(t), myR * std::sin(t), myH * t / (2.0 * M_PI));
  }
  void D1(const double t, gp_Pnt& P, gp_Vec& V) const override {
    D0(t, P);
    V.SetCoord(-myR * std::sin(t), myR * std::cos(t), myH / (2.0 * M_PI));
  }
  void D2(const double t, gp_Pnt& P, gp_Vec& V1, gp_Vec& V2) const override {
    D1(t, P, V1);
    V2.SetCoord(-myR * std::cos(t), -myR * std::sin(t), 0.0);
  }

  GeomAbs_Shape Continuity() const override { return GeomAbs_CN; }
  GeomAbs_CurveType GetType() const override { return GeomAbs_OtherCurve; }
  gp_Lin Line() const override { return gp_Lin(); }
  gp_Circ Circle() const override { return gp_Circ(); }
  gp_Elips Ellipse() const override { return gp_Elips(); }
  gp_Hypr Hyperbola() const override { return gp_Hypr(); }
  gp_Parab Parabola() const override { return gp_Parab(); }
  int Degree() const override { return 3; }
  bool IsRational() const override { return false; }
  bool IsPeriodic() const override { return false; }
  double Period() const override { return 2.0 * M_PI; }
  int NbPoles() const override { return 0; }
  int NbKnots() const override { return 0; }
  Handle(Geom_BezierCurve) Bezier() const override { return {}; }
  Handle(Geom_BSplineCurve) BSpline() const override { return {}; }
  Handle(Geom_OffsetCurve) OffsetCurve() const override { return {}; }
  int NbIntervals(const GeomAbs_Shape) const override { return 1; }
  void Intervals(TColStd_Array1OfReal& T, const GeomAbs_Shape) const override {
    T(T.Lower()) = myT0; T(T.Upper()) = myT1;
  }
  Handle(Adaptor3d_Curve) Trim(double f, double l, double) const override {
    return new HelixAdaptor(myR, myH, f, l);
  }
  gp_Vec DN(const double t, const int N) const override {
    gp_Pnt P; gp_Vec V1, V2;
    if (N == 1) { D1(t,P,V1); return V1; }
    if (N == 2) { D2(t,P,V1,V2); return V2; }
    return gp_Vec();
  }
  gp_Pnt Value(const double t) const override { gp_Pnt P; D0(t,P); return P; }

private:
  double myR, myH, myT0, myT1;
};

static void Demo5_SingleCurveRebuild()
{
  // Helix: radius=1, one full turn, height=1
  Handle(HelixAdaptor) helix = new HelixAdaptor(1.0, 1.0, 0.0, 2.0 * M_PI);

  // Rebuild with C2 continuity, tolerance 1e-4, max degree 9, max 20 segments
  Approx_Curve3d rebuilder(helix,
                           /*Tol3d=*/1.0e-4,
                           /*Order=*/GeomAbs_C2,
                           /*MaxSeg=*/20,
                           /*MaxDeg=*/9);

  if (rebuilder.IsDone())
    PrintCurveInfo(rebuilder.Curve(), "Demo5 (single curve rebuild C2)");
  else
    std::cout << "Demo5: HasResult=" << rebuilder.HasResult()
              << "  MaxError=" << rebuilder.MaxError() << "\n";
  std::cout << "  MaxError=" << rebuilder.MaxError() << "\n";
}

// ===========================================================================
// Demo 6: Curve chain rebuild
// ===========================================================================
// We have 3 Bezier arc segments forming a chain. We rebuild each segment
// individually (Approx_Curve3d per segment) then report the results.
// For G1 matching at joints, one would also enforce matching end tangents;
// here we show the independent rebuild and note the joint error.
//
// Source: a piecewise-linear-ish 3D path built from three circular arcs
// on orthogonal planes, each covering 90°.
// ===========================================================================
static void Demo6_CurveChainRebuild()
{
  // Three arcs: XY-plane, YZ-plane, XZ-plane, each 90°
  struct ArcDef { gp_Ax2 ax; double r; double t0; double t1; };
  std::vector<ArcDef> arcs = {
    { gp_Ax2(gp_Pnt(0,0,0), gp_Dir(0,0,1)), 1.0, 0.0,    M_PI/2.0 },  // XY arc
    { gp_Ax2(gp_Pnt(1,1,0), gp_Dir(1,0,0)), 1.0, -M_PI/2.0, 0.0  },  // YZ arc
    { gp_Ax2(gp_Pnt(1,0,1), gp_Dir(0,1,0)), 1.0, 0.0,    M_PI/2.0 },  // XZ arc
  };

  for (int k = 0; k < (int)arcs.size(); ++k) {
    Handle(Geom_Circle) circ = new Geom_Circle(gp_Circ(arcs[k].ax, arcs[k].r));
    Handle(Geom_TrimmedCurve) arc = new Geom_TrimmedCurve(circ, arcs[k].t0, arcs[k].t1);
    GeomAdaptor_Curve adaptor(arc);
    Handle(GeomAdaptor_Curve) hAdaptor = new GeomAdaptor_Curve(arc);

    Approx_Curve3d rebuilder(hAdaptor,
                             /*Tol3d=*/1.0e-5,
                             /*Order=*/GeomAbs_C1,
                             /*MaxSeg=*/10,
                             /*MaxDeg=*/7);

    char label[64];
    std::snprintf(label, sizeof(label), "Demo6 segment %d", k + 1);
    if (rebuilder.IsDone())
      PrintCurveInfo(rebuilder.Curve(), label);
    else
      std::cout << label << ": partial, maxErr=" << rebuilder.MaxError() << "\n";
  }
}

// ===========================================================================
// Demo 7: Single curve rebuild with discontinuity parameters
// ===========================================================================
// Source: a curve that has a sharp corner (C0 only) at a known parameter.
// We tell AdvApprox_ApproxAFunction to prefer cutting at that point via
// AdvApprox_PrefAndRec, so the BSpline knot lands exactly there.
//
// Design: two smooth arcs joined at t=π with a deliberate kink.
//   arc1: (cos t, sin t, 0)        for t ∈ [0, π]
//   arc2: (cos t, -sin t + 2, 0.5) for t ∈ [π, 2π]  (different plane → kink at π)
// ===========================================================================

class KinkCurveEvaluator : public AdvApprox_EvaluatorFunction
{
public:
  // Result layout: 3 components (X,Y,Z) — one 3D sub-space
  void operator()(int*    Numpoints,
                  double* Parameter,
                  int*    DerivativeRequest,
                  double* Result,
                  int*    ErrorCode) const override
  {
    double t  = *Parameter;
    int    d  = *DerivativeRequest;
    *ErrorCode = 0;

    auto eval = [&](double tt, double* r) {
      if (tt <= M_PI) {
        r[0] = std::cos(tt);
        r[1] = std::sin(tt);
        r[2] = 0.0;
      } else {
        r[0] = std::cos(tt);
        r[1] = -std::sin(tt) + 2.0;
        r[2] = 0.5;
      }
    };

    double h = 1.0e-6;
    if (d == 0) {
      eval(t, Result);
    } else if (d == 1) {
      double r0[3], r1[3];
      eval(t - h, r0); eval(t + h, r1);
      for (int i = 0; i < 3; ++i) Result[i] = (r1[i] - r0[i]) / (2.0 * h);
    } else if (d == 2) {
      double r0[3], r1[3], r2[3];
      eval(t - h, r0); eval(t, r1); eval(t + h, r2);
      for (int i = 0; i < 3; ++i) Result[i] = (r0[i] - 2*r1[i] + r2[i]) / (h*h);
    }
    (void)Numpoints;
  }
};

static void Demo7_CurveRebuildWithDiscontinuity()
{
  // Preferred cut at the kink parameter π; recommended at midpoints
  NCollection_Array1<double> pref(1, 1);  pref(1) = M_PI;
  NCollection_Array1<double> rec (1, 3);
  rec(1) = M_PI / 2.0;
  rec(2) = M_PI;
  rec(3) = 3.0 * M_PI / 2.0;

  AdvApprox_PrefAndRec cutter(pref, rec);
  KinkCurveEvaluator   evaluator;

  // Tolerances: one 3D sub-space
  Handle(NCollection_HArray1<double>) tol1d = new NCollection_HArray1<double>(1,0);
  Handle(NCollection_HArray1<double>) tol2d = new NCollection_HArray1<double>(1,0);
  Handle(NCollection_HArray1<double>) tol3d = new NCollection_HArray1<double>(1,1);
  tol3d->SetValue(1, 1.0e-4);

  AdvApprox_ApproxAFunction approx(
    /*Num1D=*/0, /*Num2D=*/0, /*Num3D=*/1,
    tol1d, tol2d, tol3d,
    /*First=*/0.0, /*Last=*/2.0 * M_PI,
    /*Cont=*/GeomAbs_C0,   // C0 only — kink at π
    /*MaxDeg=*/9,
    /*MaxSeg=*/20,
    evaluator, cutter);

  if (approx.IsDone()) {
    std::cout << "Demo7 (kink curve rebuild):"
              << "  degree=" << approx.Degree()
              << "  nbKnots=" << approx.NbKnots()
              << "  maxErr=" << approx.MaxError(3, 1) << "\n";
    // Check that the kink parameter π appears as a knot
    auto knots = approx.Knots();
    bool foundKink = false;
    for (int i = knots->Lower(); i <= knots->Upper(); ++i) {
      if (std::fabs(knots->Value(i) - M_PI) < 1.0e-8) { foundKink = true; break; }
    }
    std::cout << "  kink knot at π found: " << (foundKink ? "YES" : "NO") << "\n";
  } else {
    std::cout << "Demo7: HasResult=" << approx.HasResult() << "\n";
  }
}

// ===========================================================================
// Demo 8: Surface rebuild — D0 + partial derivatives via AdvApp2Var
// ===========================================================================
// Source: torus-like surface  F(u,v) = ((R + r*cos v)*cos u,
//                                       (R + r*cos v)*sin u,
//                                        r*sin v)
//   u ∈ [0, π/2],  v ∈ [0, π/2]   (one quarter patch)
//   R = 2, r = 0.5
//
// AdvApp2Var_ApproxAFunc2Var internally samples D0, and first/second
// derivatives in U and V (Duu, Duv, Dvv, surface normals) as needed
// to enforce the requested continuity.
// ===========================================================================

class TorusEvaluator : public AdvApp2Var_EvaluatorFunc2Var
{
public:
  // AdvApp2Var_EvaluatorFunc2Var::operator() signature:
  //   (OrderU, OrderV, U, V, Result, ErrorCode)
  // OrderU/V = 0 → value; 1 → first deriv; 2 → second deriv
  void operator()(int*    OrderU,
                  int*    OrderV,
                  double* U,
                  double* V,
                  double* Result,
                  int*    ErrorCode) const override
  {
    double u = *U, v = *V;
    const double R = 2.0, r = 0.5;

    *ErrorCode = 0;

    // F(u,v)
    double rho = R + r * std::cos(v);
    double x = rho * std::cos(u);
    double y = rho * std::sin(u);
    double z = r   * std::sin(v);

    if (*OrderU == 0 && *OrderV == 0) {
      Result[0] = x; Result[1] = y; Result[2] = z;
    }
    // D/Du
    else if (*OrderU == 1 && *OrderV == 0) {
      Result[0] = -rho * std::sin(u);
      Result[1] =  rho * std::cos(u);
      Result[2] =  0.0;
    }
    // D/Dv
    else if (*OrderU == 0 && *OrderV == 1) {
      double drho_dv = -r * std::sin(v);
      Result[0] = drho_dv * std::cos(u);
      Result[1] = drho_dv * std::sin(u);
      Result[2] = r * std::cos(v);
    }
    // D²/Du²
    else if (*OrderU == 2 && *OrderV == 0) {
      Result[0] = -rho * std::cos(u);
      Result[1] = -rho * std::sin(u);
      Result[2] =  0.0;
    }
    // D²/Dv²
    else if (*OrderU == 0 && *OrderV == 2) {
      double d2rho = -r * std::cos(v);
      Result[0] = d2rho * std::cos(u);
      Result[1] = d2rho * std::sin(u);
      Result[2] = -r * std::sin(v);
    }
    // D²/DuDv
    else if (*OrderU == 1 && *OrderV == 1) {
      double drho_dv = -r * std::sin(v);
      Result[0] = -drho_dv * std::sin(u);
      Result[1] =  drho_dv * std::cos(u);
      Result[2] =  0.0;
    }
    else {
      *ErrorCode = 1;  // order not supported
    }
  }
};

static void Demo8_SurfaceRebuild()
{
  Handle(NCollection_HArray1<double>) tol1d = new NCollection_HArray1<double>(1,0);
  Handle(NCollection_HArray1<double>) tol2d = new NCollection_HArray1<double>(1,0);
  Handle(NCollection_HArray1<double>) tol3d = new NCollection_HArray1<double>(1,1);
  tol3d->SetValue(1, 1.0e-4);

  // Boundary tolerances (same, no special edge treatment)
  Handle(NCollection_HArray2<double>) btol1d = new NCollection_HArray2<double>(1,1,1,4);
  Handle(NCollection_HArray2<double>) btol2d = new NCollection_HArray2<double>(1,1,1,4);
  Handle(NCollection_HArray2<double>) btol3d = new NCollection_HArray2<double>(1,1,1,4);
  for (int i = 1; i <= 4; ++i) btol3d->SetValue(1, i, 1.0e-4);

  TorusEvaluator    evaluator;
  AdvApprox_DichoCutting ucutter, vcutter;

  AdvApp2Var_ApproxAFunc2Var approx(
    /*Num1D=*/0, /*Num2D=*/0, /*Num3D=*/1,
    tol1d, tol2d, tol3d,
    btol1d, btol2d, btol3d,
    /*U0=*/0.0, /*U1=*/M_PI / 2.0,
    /*V0=*/0.0, /*V1=*/M_PI / 2.0,
    /*FavorIso=*/GeomAbs_IsoV,
    /*ContU=*/GeomAbs_C1,
    /*ContV=*/GeomAbs_C1,
    /*PrecisCode=*/2,
    /*MaxDegU=*/9,
    /*MaxDegV=*/9,
    /*MaxPatch=*/16,
    evaluator, ucutter, vcutter);

  if (approx.IsDone()) {
    Handle(Geom_BSplineSurface) surf =
      Handle(Geom_BSplineSurface)::DownCast(approx.Surface(1));
    if (!surf.IsNull()) {
      std::cout << "Demo8 (torus patch rebuild):"
                << "  degU=" << surf->UDegree()
                << "  degV=" << surf->VDegree()
                << "  nPolesU=" << surf->NbUPoles()
                << "  nPolesV=" << surf->NbVPoles()
                << "\n";
    }
    std::cout << "  maxErr3D=" << approx.MaxError(3, 1) << "\n";
    std::cout << "  avgErr3D=" << approx.AverageError(3, 1) << "\n";
  } else {
    std::cout << "Demo8: HasResult=" << approx.HasResult() << "\n";
  }
}

// ===========================================================================
// main
// ===========================================================================
int main()
{
  std::cout << "=== OCCT Approximation (Rebuild) Demos ===\n\n";

  std::cout << "--- Demo 5: Single curve rebuild (D0+D1+D2, C2) ---\n";
  Demo5_SingleCurveRebuild();

  std::cout << "\n--- Demo 6: Curve chain rebuild (3 arcs, C1 each) ---\n";
  Demo6_CurveChainRebuild();

  std::cout << "\n--- Demo 7: Curve rebuild with discontinuity at π ---\n";
  Demo7_CurveRebuildWithDiscontinuity();

  std::cout << "\n--- Demo 8: Surface rebuild (torus patch, C1 x C1) ---\n";
  Demo8_SurfaceRebuild();

  std::cout << "\nDone.\n";
  return 0;
}
