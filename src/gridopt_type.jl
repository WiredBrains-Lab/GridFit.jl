"""
GridOpt types and creation functions.

The GridOpt object contains the precalculated and intermediate values used during the optimization
"""

using NearestNeighbors,NIfTI,LinearAlgebra

include("grid_types.jl")
include("math.jl")

struct Edges{T<:Real}
	coords::Vector{Vector{T}}
	tree::KDTree
end

function Edges(coords::Vector{Vector{T}}) where T<:Real
	tree = KDTree(hcat(coords...))
	Edges(coords,tree)
end

function Edges(data::AbstractArray{T,3},affine::Matrix{K}) where {T<:Real,K<:Real}
    th = maximum(data)*0.99
    edge_ii = findall(CartesianIndices(data)) do ii
		data[ii]>0. || return false
		for x=-1:1
			for y=-1:1
				for z=-1:1
					j = [ii.I[1]+x,ii.I[2]+y,ii.I[3]+z]
					checkbounds(Bool,data,j...) && data[j...] == 0. && return true
				end
			end
		end
		return false
	end
    edge_coords = Vector{Float64}[]

    for c in edge_ii
        cx = ijk_xyz(data,collect(c.I),affine)
        push!(edge_coords,cx)
    end

    Edges(edge_coords)
end

Edges(data::NIVolume) = Edges(data.raw,getaffine(data))

Edges(fname::AbstractString) = Edges(niread(fname))

"""
    GridOpt(grid::Grid,edges::Edges,grid_coords::Dict{Int,Vector{Float64}})

Create a single object to track the grid optimization. 

    `grid``: an immutable `Grid` object that describes the grid geometry.
    `edges`: an immutable `Edges` object that contains the coordinates of edge points.
    `grid_coords`: the current best-guess of the grid coordinate points
    `fixed_points`: a `Vector` of points that should be treated as fixed (known points that
        should not be optimized)

"""
mutable struct GridOpt
	grid::Grid
	edges::Edges
	grid_coords::Dict{Int,Vector{Float64}}
    fixed_points::Dict{Int,Vector{Float64}}
end

GridOpt(grid::Grid,edges::Edges,fixed_points::Dict{Int,Vector{Float64}}=Dict{Int,Vector{Float64}}()) = GridOpt(grid,edges,Dict{Int,Vector{Float64}}(),fixed_points)

coord_dist(s::GridOpt,a::Int,b::Int) = norm(s.grid_coords[a] .- s.grid_coords[b])

function init_coords!(s::GridOpt,start_coords::Dict)
	grid_nums = s.grid.layout[:]

	grid_coords = [grid_ii(s.grid,i) for i in grid_nums]
	grid_realcoords = [[0.,i[2]-1,i[1]-1] .* s.grid.spacing for i in grid_coords]
	grid_shiftcoords = [i .+ start_coord for i in grid_realcoords]
	grid_coords = Dict{Int,Vector{Float64}}()
	for i in 1:length(grid_nums)
		grid_coords[grid_nums[i]] = grid_shiftcoords[i]
	end

	s.grid_coords = grid_coords
end

function rotate!(s::GridOpt,xa::K,ya::L,za::M) where {K<:Real,L<:Real,M<:Real}
	grid_nums = sort(s.grid.layout[:])
	grid_c = hcat([s.grid_coords[c] for c in grid_nums]...)

	grid_c = rotate(grid_c,xa,ya,za)
	for i in 1:length(grid_nums)
		s.grid_coords[grid_nums[i]] = grid_c[:,i]
	end

	return s
end

function translate!(s::GridOpt,xt::K,yt::L,zt::M) where {K<:Real,L<:Real,M<:Real}
	grid_nums = sort(s.grid.layout[:])
	grid_c = hcat([s.grid_coords[c] for c in grid_nums]...)

	grid_c = translate(grid_c,xt,yt,zt)
	for i in 1:length(grid_nums)
		s.grid_coords[grid_nums[i]] = grid_c[:,i]
	end

	return s
end
