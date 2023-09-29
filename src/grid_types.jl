"""
Types describing the different orientations of grids and functions to create and
manipulate them.
"""

using LinearAlgebra

abstract type Grid end

struct RectGrid <: Grid
	layout::Matrix{Int}
	spacing::Float64
end

RectGrid(columns::Int,rows::Int;spacing=4.0,trim=[:,:]) = RectGrid(reverse(reshape(1:(columns*rows),(columns,rows))',dims=2)[trim...],spacing)

struct SquareGrid <: Grid
	layout::Matrix{Int}
	spacing::Float64
end

SquareGrid(width::Int;spacing=4.0,trim=[:,:]) = SquareGrid(reverse(reshape(1:(width^2),(width,width))',dims=2)[trim...],spacing)

struct ArbitraryGrid <: Grid
	layout::Matrix{Int}
	spacing::Float64
end

grid_ii(g::Grid,i::Int) = collect(findfirst(isequal(i),g.layout).I)
grid_dist(g::Grid,a::Int,b::Int) = norm(grid_ii(g,a) .- grid_ii(g,b)) * g.spacing