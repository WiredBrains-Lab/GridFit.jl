include("gridopt_type.jl")

using NIfTI

function make_dot(dset,ii::Vector{T},r::Int,v::Real=1) where T<:Real
    for x=-r:r
        for y=-r:r
            for z=-r:r
                if x^2+y^2+z^2<=r^2
					c = [ii[1]+x,ii[2]+y,ii[3]+z]
					if dset[c...]==0
						dset[c...] = v
					else
						dset[c...] = mean([v,dset[c...]])
					end
                end
            end
        end
    end
end

function writegrid(out_dset::NIVolume,out_dat,grid_coords::Dict,fname::String,values::Union{Real,Dict}=1.;r::Int=2)
	for (k,c) in grid_coords
		if isa(values,Real)
	    	make_dot(out_dat,xyz_ijk(out_dset,c),r,values)
		else
			haskey(values,k) && make_dot(out_dat,xyz_ijk(out_dset,c),r,values[k])
		end
	end

	niwrite(datadir(fname),NIVolume(out_dset.header,out_dset.extensions,out_dat))
end

function writegrid(anat_fname::String,grid_coords::Dict,fname::String,values::Union{Real,Dict}=1.;r::Int=2)
	if isfile(fname)
		writegrid(grid_coords,fname,values;r)
	else
		out_dset = niread(anat_fname)
		out_dat = zeros(Float32,size(out_dset))

		writegrid(out_dset,out_dat,grid_coords,fname,values;r)
	end
end

"""
	function writegrid(grid_coords::Dict,fname::String,values::Union{Real,Dict}=1.;r::Float64=2.)

Writes the given coordinates as spheres of radius `r` and value `values` to a NIfTI file `fname`. Currently not reliable,
and I'm not sure why. I don't know enough about the NIfTI file format to fully debug this. For now I'm using the alternative
function `writegrid_afni`, which just passes to coordinates to an AFNI program.
"""
function writegrid(grid_coords::Dict,fname::String,values::Union{Real,Dict}=1.;r::Float64=2.)
	out_dset = niread(fname)
	out_dat = Float32.(out_dset.raw)

	writegrid(out_dset,out_dat,grid_coords,fname,values;r)
end

"""
	function writegrid_afni(anat_fname::String,grid_coords::Dict,fname::String,values::Union{Real,Dict}=1.;r::Float64=2.)

Similar to `writegrid`, but uses the AFNI program `3dUndump` instead of the NIfTI.jl library. I've found this function
to be more reliable (although it has an external dependency).
"""
function writegrid_afni(anat_fname::String,grid_coords::Dict,fname::String,values::Union{Real,Dict}=1.;r::Float64=2.)
	vals = []
	for (k,c) in grid_coords
		if isa(values,Real)
			push!(vals,c)
		else
			haskey(values,k) && push!(vals,vcat(c,values[k]))
		end
	end

	dval = isa(values,Real) ? values : 1

	cmd = `3dUndump -prefix $fname -master $anat_fname -dval $dval -xyz -srad $r -orient LPI -`

	open(pipeline(cmd),"w",stdout) do f
		for v in vals
			write(f,join(v," "))
			write(f,"\n")
		end
	end
end