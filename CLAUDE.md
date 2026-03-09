
please write markdown document for the code about the fitting.

code in the 
```txt
    src/ModelingData/TKG3d
        AdvApprox       // 1d curve fitting with adaptive segment, including eval sampling, cutting, ..., 
                        // and single segment fitting AdvApprox_SimpleApprox
                        // please notes,  single segment use both jacobi approx and hermite interp.
                        // because inter-segment joint need maintain
    src/ModelingData/TKGeomBase
        AdvApp2Var      // 2d surf fitting based on 1d curve fitting
        AppCont         // ?
        AppDef          // do fitting using AdvApprox/AdvApp2Var or by self, usually with direct solve or gradient/BFGS method, some method add regularization term to make curve smooth.
                        // the data is AppDef_MultiPointConstraint: which is point from multiple curve at same u.  AppDef_MultiLine is multiple line data, each u with AppDef_MultiPointConstraint
        AppParCurves    // ? seems to base class of AppDef
        Approx          // final wrapper

    src/ModelingAlgorithms
        TKTopAlgo       // brep topology fitting
```

    now write markdown document to mydoc/
    each directory one markdown file.

    write:
    1. architecture
    2. purpose of the module
    3. simple principle, theory, algorithm of the module
    4. explain file by file, purpose of each file, and simple implement, key algorithm or math theory
    5. finally give summary table of the files
   

    then design test solution to write demo example to mydemo:

    interp:
    1. curve interp with points
    2. curve interp with points and end tangent
    3. curve interp with points and tangent
    4. surf interp with points, uv, normal
    
    approx:
    5. single curve rebuild, which may sample D0 position, D1 tangent, D2 curvature
    6. curve chain rebuild, which may sample D0 position, D1 tangent, D2 curvature
    7. single curve rebuild, which may sample D0 position, D1 tangent, D2 curvature, with discontinuity points array
    8. surf rebuild, which may sample D0 position, D1 tangent, D2 curvature, Duu Duv Dvv normal.

    the demo curve and surf, you should design it.


