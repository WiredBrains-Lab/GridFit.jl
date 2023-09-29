"""
Simple math functions.
"""

using LinearAlgebra,NIfTI

function angle(a::AbstractVector{T},b::AbstractVector{K},c::AbstractVector{L}) where {T<:Real,K<:Real,L<:Real}
	@assert length(a) == length(b) == length(c)
	ab = a .- b
	cb = c .- b
	acos(dot(normalize(ab),normalize(cb)))
end

function rotate(coords::Matrix{T},xa::K,ya::L,za::M) where {T<:Real,K<:Real,L<:Real,M<:Real}
	Rx = [1. 0. 0.; 0. cos(xa) -sin(xa); 0. sin(xa) cos(xa)]
	Ry = [cos(ya) 0. sin(ya); 0. 1. 0.; -sin(ya) 0. cos(ya)]
	Rz = [cos(za) -sin(za) 0.; sin(za) cos(za) 0.; 0. 0. 1.]

	c = Rx*coords
	c = Ry*c
	c = Rz*c

	return c
end

function translate(coords::Matrix{T},xt::K,yt::L,zt::M) where {T<:Real,K<:Real,L<:Real,M<:Real}
	c = copy(coords)
	c[1,:] .+= xt
	c[2,:] .+= yt
	c[3,:] .+= zt
	return c
end

function ijk_xyz(dset::AbstractArray{T,3},ijk::Vector{Int},affine::Matrix{K}) where {T<:Real,K<:Real}
    @assert length(ijk)==3
    size(affine)==(3,4) && (d = vcat(d,[0. 0. 0. 1.]))
    @assert size(affine)==(4,4)

    return (affine * vcat(ijk .- 1,[1]))[1:3]
end

ijk_xyz(dset::NIVolume,ijk::Vector{Int}) = ijk_xyz(dset.raw,ijk,getaffine(dset))