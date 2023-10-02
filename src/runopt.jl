"""
Functions to perform the actual grid optimization
"""

include("gridopt_type.jl")

using Optim

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

function fixed_dist(s::GridOpt)
	dists = Float64[]

    for (p,v) in s.fixed_points
		if haskey(s.grid_coords,p)
			push!(dists,norm(s.grid_coords[p] .- v))
		end
	end
	mean(dists)
end


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
function rigid_dist(s::GridOpt,x::Vector{Float64},a)
	ss = deepcopy(s)
	translate!(ss,x[1:3]...)
	rotate!(ss,x[4:6]...)

	push!(a,copy(ss.grid_coords))
	return fixed_dist(ss)
end

"""
    full_dist(...)

A combination of several cost functions to optimize several at once.
"""
function full_dist(s::GridOpt,i::Int,x::Vector{Float64};neighbor=15.,fixed=10.,edge=1.)
	ss = deepcopy(s)

	ss.grid_coords[i] .+= x

	return neighbor*neighbor_dist(ss,i) + fixed*fixed_dist(ss) + edge*edge_dist(ss,i) + angle_dist(ss,i)
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
    rigid_fit!(s::GridOpt)

Performs a rigid fit (doesn't move any electrodes relative to
one another), based on the cost of the fixed points.

This function works well as a first approximation to get the 
grid in the approximate location based on the fixed points.
"""
function rigid_fit!(s::GridOpt)
	saves = []
	o = optimize(x->rigid_dist(s,x,saves),zeros(6),Optim.Options(iterations=3_000_000,store_trace=true))
	x = Optim.minimizer(o)
	translate!(s,x[1:3]...)
	rotate!(s,x[4:6]...)
	return saves
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
    run_gridopt!(s::GridOpt;iters::Int=20_000,change_thresh=1e-6,change_count=20,verbose=true)

Run the grid optimization on the `GridOpt` object.
"""
function run_gridopt!(s::GridOpt;iters::Int=20_000,change_thresh=1e-6,change_count=20,verbose=true)
	final_coords = Dict()

    verbose && @info "Starting grid optimization"
    last_d = 1e6
    change_c = 0
    verbose && @info "- Rigid fitting starting coordinates"
    o = rigid_fit!(s)
    verbose && @info "- Optimizing individual electrodes:"
    for i=1:iters
        verbose && @info "    - iter $i"
        iterate!(s)
        dd = full_dist(s)
        verbose && @info dd
        perc_change = 100*(last_d - dd.total) / last_d
        verbose && @info perc_change
        abs(perc_change) < change_thresh ? (change_c += 1) : (change_c = 0)
        change_c>change_count && break
        last_d = dd.total
    end

    return s
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
    init_coords!(s,first(values(fixed_points)))

    run_gridopt!(s;args...)
end