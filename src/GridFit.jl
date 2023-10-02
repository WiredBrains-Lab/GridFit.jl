module GridFit

export GridOpt,Edges
export RectGrid,SquareGrid,ArbitraryGrid
export run_gridopt,run_gridopt!
export writegrid_afni

include("runopt.jl")
include("output.jl")

end
