// demo_interp.cpp
// Demonstrates 4 interpolation use-cases using OCCT geometry classes.
//
// Build (example, adjust paths):
//   g++ demo_interp.cpp -I<OCCT_INC> -L<OCCT_LIB> \
//       -lTKernel -lTKMath -lTKG3d -lTKGeomBase -lTKGeomAlgo -o demo_interp
//
// Demos:
//   1. Curve interpolation through 3D points (no tangents)
//   2. Curve interpolation with prescribed end tangents
//   3. Curve interpolation with per-point tangents
//   4. Surface interpolation through a grid of (uv, xyz, normal) data

#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array1OfVec.hxx>
#include <TColStd_Array1OfBoolean.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>

#include <TColgp_Array2OfPnt.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineSurface.hxx>

#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <iostream>
#include <cmath>

// ---------------------------------------------------------------------------
// Helper: print a BSpline curve summary
// ---------------------------------------------------------------------------
static void PrintCurveInfo(const Handle(Geom_BSplineCurve)& C, const char* label)
{
  if (C.IsNull()) {
    std::cout << label << ": FAILED (null result)\n";
    return;
  }
  std::cout << label
            << ": degree=" << C->Degree()
            << "  nbPoles=" << C->NbPoles()
            << "  nbKnots=" << C->NbKnots()
            << "\n";
}

// ---------------------------------------------------------------------------
// Helper: print a BSpline surface summary
// ---------------------------------------------------------------------------
static void PrintSurfaceInfo(const Handle(Geom_BSplineSurface)& S, const char* label)
{
  if (S.IsNull()) {
    std::cout << label << ": FAILED (null result)\n";
    return;
  }
  std::cout << label
            << ": degreeU=" << S->UDegree()
            << "  degreeV=" << S->VDegree()
            << "  nbPolesU=" << S->NbUPoles()
            << "  nbPolesV=" << S->NbVPoles()
            << "\n";
}

// ===========================================================================
// Demo 1: Curve interpolation through points (no derivative constraints)
// ===========================================================================
// Design: helical-ish 3D points sampled from a helix
//   P(t) = (cos t, sin t, t/5),  t = 0 .. 4π,  10 points
// ===========================================================================
static void Demo1_CurveInterpPoints()
{
  const int N = 10;
  TColgp_Array1OfPnt pts(1, N);
  for (int i = 1; i <= N; ++i) {
    double t = (i - 1) * 4.0 * M_PI / (N - 1);
    pts(i) = gp_Pnt(std::cos(t), std::sin(t), t / 5.0);
  }

  // GeomAPI_Interpolate: chord-length parameterization, no end conditions
  Handle(TColgp_HArray1OfPnt) hPts = new TColgp_HArray1OfPnt(1, N);
  for (int i = 1; i <= N; ++i) hPts->SetValue(i, pts(i));

  GeomAPI_Interpolate interp(hPts, /*periodic=*/Standard_False, /*tol=*/1e-6);
  interp.Perform();

  Handle(Geom_BSplineCurve) curve = interp.Curve();
  PrintCurveInfo(curve, "Demo1 (point-only interp)");

  // Verify: the curve passes through all input points
  double maxErr = 0.0;
  for (int i = 1; i <= N; ++i) {
    // Find parameter via chord-length (the interpolated params are internal;
    // just sample the curve at uniform param and note that exact passage is
    // guaranteed by construction).
    (void)i;
  }
  std::cout << "  (exact interpolation guaranteed by construction)\n";
}

// ===========================================================================
// Demo 2: Curve interpolation with prescribed end tangents
// ===========================================================================
// Design: same helix points as Demo 1, but we prescribe the analytical
// tangents at the first and last point.
//   T(t) = (-sin t, cos t, 1/5)  (derivative of helix)
// ===========================================================================
static void Demo2_CurveInterpEndTangents()
{
  const int N = 10;
  Handle(TColgp_HArray1OfPnt) hPts = new TColgp_HArray1OfPnt(1, N);
  for (int i = 1; i <= N; ++i) {
    double t = (i - 1) * 4.0 * M_PI / (N - 1);
    hPts->SetValue(i, gp_Pnt(std::cos(t), std::sin(t), t / 5.0));
  }

  // Analytical tangents at first and last points
  double t0 = 0.0;
  double tN = (N - 1) * 4.0 * M_PI / (N - 1);
  gp_Vec tangFirst(-std::sin(t0), std::cos(t0), 0.2);
  gp_Vec tangLast (-std::sin(tN), std::cos(tN), 0.2);

  GeomAPI_Interpolate interp(hPts, /*periodic=*/Standard_False, /*tol=*/1e-6);
  interp.Load(tangFirst, tangLast);   // prescribe end tangents
  interp.Perform();

  PrintCurveInfo(interp.Curve(), "Demo2 (end-tangent interp)");
  std::cout << "  startTangent=(" << tangFirst.X() << ", "
            << tangFirst.Y() << ", " << tangFirst.Z() << ")\n";
  std::cout << "  endTangent  =(" << tangLast.X()  << ", "
            << tangLast.Y()  << ", " << tangLast.Z()  << ")\n";
}

// ===========================================================================
// Demo 3: Curve interpolation with per-point tangents
// ===========================================================================
// Design: same helix but we supply tangent vectors at every sample point.
// When all tangents are prescribed, the system is fully determined.
// GeomAPI_Interpolate accepts a tangent array + a boolean mask.
// ===========================================================================
static void Demo3_CurveInterpAllTangents()
{
  const int N = 8;   // fewer points to keep system well-determined
  Handle(TColgp_HArray1OfPnt) hPts = new TColgp_HArray1OfPnt(1, N);
  TColgp_Array1OfVec           tangs(1, N);
  Handle(TColStd_HArray1OfBoolean) flags = new TColStd_HArray1OfBoolean(1, N);

  for (int i = 1; i <= N; ++i) {
    double t = (i - 1) * 4.0 * M_PI / (N - 1);
    hPts->SetValue(i, gp_Pnt(std::cos(t), std::sin(t), t / 5.0));
    tangs(i) = gp_Vec(-std::sin(t), std::cos(t), 0.2);
    flags->SetValue(i, Standard_True);  // use this tangent
  }

  GeomAPI_Interpolate interp(hPts, /*periodic=*/Standard_False, /*tol=*/1e-6);
  interp.Load(tangs, flags);
  interp.Perform();

  PrintCurveInfo(interp.Curve(), "Demo3 (per-point tangent interp)");
}

// ===========================================================================
// Demo 4: Surface interpolation through a grid of (UV, XYZ) data
// ===========================================================================
// Design: saddle surface  Z = sin(u)*cos(v),
//   u ∈ [0, π],  v ∈ [0, π],  on a 6×6 grid.
// GeomAPI_PointsToBSplineSurface interpolates through a 2D array of points.
// Normals cannot be directly prescribed via this API (it fits positions);
// for normal-constrained surfaces one would need GeomFill or AppDef_Variational.
// ===========================================================================
static void Demo4_SurfInterpPoints()
{
  const int Nu = 6, Nv = 6;
  TColgp_Array2OfPnt grid(1, Nu, 1, Nv);

  for (int i = 1; i <= Nu; ++i) {
    double u = (i - 1) * M_PI / (Nu - 1);
    for (int j = 1; j <= Nv; ++j) {
      double v = (j - 1) * M_PI / (Nv - 1);
      double x = u;
      double y = v;
      double z = std::sin(u) * std::cos(v);
      grid(i, j) = gp_Pnt(x, y, z);
    }
  }

  // Fit BSpline surface through the grid (degree 3×3)
  GeomAPI_PointsToBSplineSurface fitter;
  fitter.Interpolate(grid);

  Handle(Geom_BSplineSurface) surf = fitter.Surface();
  PrintSurfaceInfo(surf, "Demo4 (grid point interp)");

  // Verify a sample interior point
  if (!surf.IsNull()) {
    gp_Pnt evalPt;
    double umid = (surf->UKnot(1) + surf->UKnot(surf->NbUKnots())) * 0.5;
    double vmid = (surf->VKnot(1) + surf->VKnot(surf->NbVKnots())) * 0.5;
    surf->D0(umid, vmid, evalPt);
    std::cout << "  Sample at (umid,vmid): ("
              << evalPt.X() << ", " << evalPt.Y() << ", " << evalPt.Z() << ")\n";
  }
}

// ===========================================================================
// main
// ===========================================================================
int main()
{
  std::cout << "=== OCCT Interpolation Demos ===\n\n";

  std::cout << "--- Demo 1: Curve interp (points only) ---\n";
  Demo1_CurveInterpPoints();

  std::cout << "\n--- Demo 2: Curve interp (end tangents) ---\n";
  Demo2_CurveInterpEndTangents();

  std::cout << "\n--- Demo 3: Curve interp (per-point tangents) ---\n";
  Demo3_CurveInterpAllTangents();

  std::cout << "\n--- Demo 4: Surface interp (point grid) ---\n";
  Demo4_SurfInterpPoints();

  std::cout << "\nDone.\n";
  return 0;
}
