# GridFit

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://azraq27.github.io/GridFit.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://azraq27.github.io/GridFit.jl/dev/)
[![Build Status](https://github.com/azraq27/GridFit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/azraq27/GridFit.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/azraq27/GridFit.jl.svg?branch=main)](https://travis-ci.com/azraq27/GridFit.jl)
[![Coverage](https://codecov.io/gh/azraq27/GridFit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/azraq27/GridFit.jl)

This package is designed to extrapolate coordinates from an electrode grid placed over the surface of the brain (a common problem in ECoG neuroscience research). It requires a skull-stripped MRI (to define the pial surface), the electrode grid geometry, and the coordinates of known points (in the coordinate system of the MRI).

I initially tried to find an analytic geometric solution to this, but then decided to just stay simple and solve the problem with CPU time. This algorithm uses a simple cost function based on the known factors and a simple non-linear optimization algorithm to find the best electrode coordinates.

## Example:

    using GridFit

Filename to the already skullstripped MRI:

    anat_dset = "anatomy_synthstrip.nii.gz"

This grid has 64-contacts in a 16x4 arrangement with 4mm spacing:

    grid = RectGrid(16,4,spacing=4.0)

If you're wondering how the electrodes are numbered in this layout, print out the `layout` property:

    julia> grid.layout
    4×16 Matrix{Int64}:
     16  15  14  13  12  11  10   9   8   7   6   5   4   3   2   1
     32  31  30  29  28  27  26  25  24  23  22  21  20  19  18  17
     48  47  46  45  44  43  42  41  40  39  38  37  36  35  34  33
     64  63  62  61  60  59  58  57  56  55  54  53  52  51  50  49

These are the measure coordinates (from the intraoperative navigation system on the exposed electrodes, in the coordinate system of the anatomy MRI:

    coords = Dict(
         16 => [-67.83, -0.546, -24.928],     
         32 => [-66.837, -2.577, -27.328],
         48 => [-67.853, -3.085, -33.728],
         64 => [-67.853, -4.609, -37.728],
         63 => [-67.853, -9.179, -36.928],
         62 => [-69.376, -14.257, -36.128],
         61 => [-67.853, -18.827, -33.728],
         60 => [-66.837, -21.874, -32.928]
    )

Run the optimization:

    opt = run_gridopt(grid,anat_dset,coords)

The `run_gridopt` function first translate/rotates the entire grid to be close to the fixed coordinates ("rigid fitting") and then iteratively goes through each electrode and tries to optimize its location simultaneously based on:
- `neighbor`: distance to its neighbors (using the grid spacing)
- `edge`: distance to the pial surface
- `fixed`: distance to the known coordinates (only if this is one of the given fixed coordinates)
- `angle`: angle with surrounding electrodes (using the grid geometry)

On each iteration it will print the total loss, percentage change, and continue optimizing until the loss stops decreasing.

    [ Info: Starting grid optimization
    [ Info: - Rigid fitting starting coordinates
    [ Info: - Optimizing individual electrodes:
    [ Info:     - iter 1
    [ Info: (neighbor = 1.000149305829309, edge = 2.9928794459057793, fixed = 1.5574529108189121, angle = 1.0000473415185291, total = 34.56969548305307)
    [ Info: 99.99654303045169
    [ Info:     - iter 2
    [ Info: (neighbor = 1.0002292302432017, edge = 2.9764837140212506, fixed = 1.5502882409997443, angle = 1.0000808335078437, total = 34.482885411174564)
    [ Info: 0.2511161023129621
    [ Info:     - iter 3
    [ Info: (neighbor = 1.0002520960698154, edge = 2.960812986274054, fixed = 1.5459388954463982, angle = 1.0000846125289367, total = 34.42406799431421)
    [ Info: 0.17056988172253357
    ...
    [ Info:     - iter 292                                                          
    [ Info: (neighbor = 1.0119593444413224, edge = 0.30534186860100215, fixed = 1.5833514303514304, angle = 1.003742645161798, total = 32.32198898389694)
    [ Info: 0.0
    [ Info:     - iter 293                                                          
    [ Info: (neighbor = 1.0119593444413224, edge = 0.30534186860100215, fixed = 1.5833514303514304, angle = 1.003742645161798, total = 32.32198898389694)
    [ Info: 0.0

The optimization will return a `GridOpt` object containing the final coordinates.

    julia> opt.grid_coords
    Dict{Int64, Vector{Float64}} with 64 entries:
      5  => [-63.1177, -44.2725, -14.2401]
      56 => [-64.2655, -35.7163, -28.8523]
      55 => [-63.7511, -39.546, -27.8198]
      35 => [-58.1649, -52.7779, -19.8579]
      60 => [-66.5365, -20.4951, -33.1204]
      30 => [-67.7287, -10.6936, -27.5931]
      32 => [-67.7698, -2.98924, -29.7472]
      6  => [-64.7875, -40.6602, -15.3284]
      45 => [-67.3333, -15.615, -30.3647]
      64 => [-67.531, -5.14021, -37.448]
      13 => [-67.3625, -13.4582, -22.6677]
      4  => [-61.2305, -47.7118, -13.1504]
      63 => [-67.5424, -8.99259, -36.3719]
      54 => [-62.7946, -43.2748, -26.7803]
      62 => [-67.8742, -12.8498, -35.2938]
      58 => [-64.8765, -28.0121, -30.9203]
      52 => [-59.0375, -50.0103, -24.7155]
      12 => [-67.1584, -17.3065, -21.5968]
      28 => [-67.248, -18.4019, -25.4431]
      ⋮  => ⋮
