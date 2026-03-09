# Geometry-Only Project — What to Take vs Skip

## Take — pure geometry, no BRep dependency

| Module | Location |
|---|---|
| `AdvApprox` | `src/ModelingData/TKG3d/AdvApprox/` |
| `AdvApp2Var` | `src/ModelingData/TKGeomBase/AdvApp2Var/` |
| `AppParCurves` | `src/ModelingData/TKGeomBase/AppParCurves/` |
| `AppCont` | `src/ModelingData/TKGeomBase/AppCont/` |
| `AppDef` | `src/ModelingData/TKGeomBase/AppDef/` |
| `Approx` | `src/ModelingData/TKGeomBase/Approx/` |

## Skip — BRep topology dependency

| Module | Why skip |
|---|---|
| `BRepApprox` | `src/ModelingAlgorithms/TKTopAlgo/BRepApprox/` — depends on BRep faces, topology, intersection results |

## Also Narrow by Use Case

If you only need **curve/surface rebuild** from a continuous function:
- Take only `AdvApprox` + `Approx_Curve3d` / `Approx_SweepApproximation`
- Skip `AppDef` entirely (that is for discrete points)

If you need **point fitting** (discrete data → BSpline):
- Take `AppParCurves` + `AppDef`
- Skip `AdvApprox` (that is for continuous functions)

If you need **both**, take everything except `BRepApprox`.
