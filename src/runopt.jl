"""
Functions to perform the actual grid optimization
"""

include("gridopt_type.jl")

using Optim,Statistics

# Cost functions:

function neighbor_dist(s::GridOpt,i::Int,r::Int=2)
	dists = Float64[]

	for x = -r:r
        for y = -r:r
			x==0 && y==0 && continue
			neighbor_ii = grid_ii(s.grid,i) .+ [y,x]
			if checkbounds(Bool,s.grid.layout,neighbor_ii...)
				neighbor_i = s.grid.layout[neighbor_ii...]
				haskey(s.grid_coords,neighbor_i) || continue

				push!(dists,norm(coord_dist(s,i,neighbor_i) .- grid_dist(s.grid,i,neighbor_i)))
			end
		end
    end

	mean(2 .^ dists)
end

function fixed_dist(fixed_points::Dict,grid_coords::Dict)
	dists = Float64[]

    for (p,v) in fixed_points
		if haskey(grid_coords,p)
			push!(dists,norm(grid_coords[p] .- v))
		end
	end
	mean(dists)
end

fixed_dist(s::GridOpt) = fixed_dist(s.fixed_points,s.grid_coords)

function edge_dist(s::GridOpt,i::Int)
	(ii,dd) = knn(s.edges.tree,s.grid_coords[i],1)
	return only(dd)
end


function angle_dist(s::GridOpt,i::Int)
	dists = Float64[]
	g = s.grid.layout
	for i=1:size(g,1)
		for j=1:size(g,2)
			b = s.grid_coords[g[i,j]]
			for neighbor in [
				([0 , -1],[-1, 0]),
				([-1,  0],[0, 1]),
				([0 ,  1],[1, 0]),
				([ 1,  0],[0, -1])
				]
				a_ii = i + neighbor[1][1]
				a_jj = j + neighbor[1][2]
				c_ii = i + neighbor[2][1]
				c_jj = j + neighbor[2][2]
				(a_ii<1 || a_ii>size(g,1)) && continue
				(c_ii<1 || c_ii>size(g,1)) && continue
				(a_jj<1 || a_jj>size(g,2)) && continue
				(c_jj<1 || c_jj>size(g,2)) && continue

				a = s.grid_coords[g[a_ii,a_jj]]
				c = s.grid_coords[g[c_ii,c_jj]]

				d = angle(a,b,c)
				push!(dists,norm(d - π/2))
			end
		end
	end
	return mean(2 .^ dists)
end

"""
    rigid_dist(...)

Calculate `fixed_dist` after applying a rigid transformation.
"""
function rigid_dist(s::GridOpt,x::Vector{Float64},a::Vector{Float64})
	rot_coords = rotate(s,x[1:3]...,a)

	return fixed_dist(s.fixed_points,rot_coords)
end

"""
    full_dist(...)

A combination of several cost functions to optimize several at once.
"""
function full_dist(s::GridOpt,i::Int,x::Vector{Float64};neighbor=15.,fixed=10.,edge=1.)
	s.grid_coords[i] .+= x

	dist = neighbor*neighbor_dist(s,i) + fixed*fixed_dist(s) + edge*edge_dist(s,i) + angle_dist(s,i)

    s.grid_coords[i] .-= x
    return dist
end

function full_dist(s::GridOpt;neighbor=15.,fixed=10.,edge=1.)
	fixed_d = fixed_dist(s)
	dd = map(collect(keys(s.grid_coords))) do i
		(neighbor=neighbor_dist(s,i), edge=edge_dist(s,i),angle=angle_dist(s,i))
	end

	neighbor_d=mean(getproperty.(dd,:neighbor))
	edge_d=mean(getproperty.(dd,:edge))
	angle=mean(getproperty.(dd,:angle))
	
	return (;neighbor=neighbor_d,edge=edge_d,fixed=fixed_d,angle,
			total=neighbor*neighbor_d + fixed*fixed_d + edge*edge_d + angle)
end

"""
	dist_limits(dist;limits...)

Returns `true` if the cost functions are in a range outside of normal values. The
default normal values were derived by running a large number of iterations and finding
the 95% bounds on each.
"""
function dist_limits(dist;
		edge = 5.577,
		angle = 1.009,
		fixed = 12.668,
		neighbor = 1.064
	)
	if  dist.edge>=edge ||
		dist.angle>=angle ||
		dist.fixed>=fixed ||
		dist.neighbor>=neighbor
		return true
	end
	return false
end

"""
    rigid_fit!(s::GridOpt,fixed::Vector{Float64})

Performs a rigid fit (doesn't move any electrodes relative to
one another), by rotating around the given start point.

This function works well as a first approximation to get the 
grid in the approximate location based on the given fixed point.
"""
function rigid_fit!(s::GridOpt)
	o = optimize(x->rigid_dist(s,x,s.start_coord),zeros(3),GradientDescent())
	x = Optim.minimizer(o)
	rotate!(s,x[1:3]...,s.start_coord)
    return nothing
end


function iterate!(s::GridOpt;neighbor=1.,fixed=1.,edge=1.)
	ii = collect(keys(s.grid_coords))
	res = Channel(length(ii))
	Threads.@threads for i in ii
		o = optimize(x->full_dist(s,i,x),zeros(3),Optim.Options(iterations=10))
		put!(res,(i,Optim.minimizer(o)))
	end
	while(isready(res))
		(i,o) = take!(res)
		s.grid_coords[i] .+= o
	end
    return nothing
end

"""
    run_gridopt!(s::GridOpt;iters=2_000,attempts=20,change_thresh=1e-6,change_count=20,verbose=true)

Run the grid optimization on the `GridOpt` object.

This function will try `iters` iterations until the cost function stops decreasing by a percent change
less than `change_thresh`. If it fails, it will try again for `attempts` times.
"""
function run_gridopt!(s::GridOpt;iters::Int=2_000,attempts::Int=20,change_thresh=1e-6,change_count=20,verbose=true)
	final_coords = Dict()

	for attempt=1:attempts
		verbose && @info "Starting grid optimization (attempt $attempt)"
		last_d = 1e6
		change_c = 0
		cost_fails = 0
		succeeded = false
		verbose && @info "- Rigid fitting starting coordinates"
		o = rigid_fit!(s)
		verbose && @info "- Optimizing individual electrodes:"
		for i=1:iters
			verbose && @info "    - iter $i"
			iterate!(s)
			dd = full_dist(s)
			verbose && @info dd
			dist_limits(dd) && (cost_fails += 1)
			if cost_fails > 20
				verbose && @info "    *** Failing due to abnormally high cost functions"
				break
			end
			perc_change = 100*(last_d - dd.total) / last_d
			verbose && @info perc_change
			abs(perc_change) < change_thresh ? (change_c += 1) : (change_c = 0)
			if change_c>change_count
				verbose && @info "Success!!"
				succeeded = true
				break
			end
			last_d = dd.total
		end

		succeeded && break
	end
	if succeeded == false
		verbose && @info "Failed to achieve successful fit. Please try again, or tweak the parameters."
		return nothing
	else
		return s
	end
end

function run_gridopt(s::GridOpt;args...)
    ss = deepcopy(s)
    run_gridopt!(ss;args...)
end

"""
    run_gridopt!(g::Grid,edge_dset::String,fixed_points::Dict;iters::Int=20_000,change_thresh=1e-6,change_count=20,verbose=true)

Run the grid optimization with the `Grid`, the `Edges` derived from the dataset `edge_dset`, 
and the known coordinates contained in `fixed_points`, a `Dict` that associates electrode 
numbers with a 3-item coordinate Vector.
"""
function run_gridopt(g::Grid,edge_dset::String,fixed_points::Dict{Int,Vector{Float64}};args...)
    edges = Edges(edge_dset)
    s = GridOpt(g,edges,fixed_points)
    init_coords!(s)

    run_gridopt!(s;args...)
end
