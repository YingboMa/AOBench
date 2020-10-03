# This file is translated from `aobench.c` by Syoyo Fujita
module AOBench

export Vect, main

using LinearAlgebra, UnPack, Setfield, MuladdMacro, Images, FileIO
using Base.Cartesian

struct Vect{T} <: AbstractVector{T}
    x::T
    y::T
    z::T
end
Vect(x, y, z) = Vect(promote(x, y, z)...,)
Base.getindex(v::Vect, i::Int) = getfield(v, i)
Base.size(::Vect) = (3,)
Base.Tuple(v::Vect) = v.x, v.y, v.z
Vect(xx::Tuple) = Vect(xx...,)
Vect(x::Number) = Vect(@ntuple 3 i->x)

struct Intersect{T,V}
    t::T
    p::V
    n::V
    hit::Int
end

struct Sphere{V,T}
    center::V
    radius::T
end

struct Plane{V}
    p::V
    n::V
end

struct Ray{V}
    org::V
    dir::V
end

struct State{T,V}
    spheres::NTuple{3,Sphere{V,T}}
    plane::Plane{V}
    width      ::Int
    height     ::Int
    nsubsamples::Int
    nao_samples::Int
end

@muladd LinearAlgebra.dot(v0::Vect, v1::Vect) = v0.x * v1.x + v0.y * v1.y + v0.z * v1.z
@muladd LinearAlgebra.cross(v0::Vect, v1::Vect) = Vect(
    v0.y * v1.z - v0.z * v1.y,
    v0.z * v1.x - v0.x * v1.z,
    v0.x * v1.y - v0.y * v1.x,
   )
LinearAlgebra.norm(v::Vect) = sqrt(dot(v, v))
LinearAlgebra.normalize(v::Vect) = v * inv(norm(v))
for op in [:+, :-, :*, :/]
    @eval Base.$op(v0::Vect, v1::Vect) = Vect(@ntuple 3 i -> $op(v0[i], v1[i]))
    @eval Base.$op(v0::Vect, s::Number) = Vect(@ntuple 3 i -> $op(v0[i], s))
    if op !== :/
        @eval Base.$op(s::Number, v0::Vect) = Vect(@ntuple 3 i -> $op(s, v0[i]))
    end
end
Base.:(-)(v::Vect{T}) where T = Vect(@ntuple 3 i -> -v[i])
Base.muladd(v1::Vect, v2::Vect, v3::Vect) = Vect(@ntuple 3 i -> muladd(v1[i], v2[i], v3[i]))
Base.zero(::Type{Vect{T}}) where T = Vect(@ntuple 3 _->zero(T))
Base.zero(::T) where {T<:Vect} = zero(T)

@fastmath @muladd function ray_intersect(isect::Intersect, ray::Ray, sphere::Sphere)
    rs = ray.org - sphere.center
    B = rs'ray.dir
    C = rs'rs - sphere.radius^2
    D = B^2 - C
    if D > zero(D)
        t = -B - sqrt(D)
        if t > zero(t) && t < isect.t
            @set! isect.t = t
            @set! isect.hit = 1

            @set! isect.p = ray.org + ray.dir * t
            @set! isect.n = normalize(isect.p - sphere.center)
        end
    end
    return isect
end

@muladd function ray_intersect(isect::Intersect, ray::Ray, plane::Plane)
    d = -plane.p'plane.n
    v = ray.dir'plane.n

    abs(v) < 1.0e-17 && return isect

    t = -(ray.org'plane.n + d) / v

    if t > zero(t) && t < isect.t
        @set! isect.t = t
        @set! isect.hit = 1

        @set! isect.p = ray.org + ray.dir * t
        @set! isect.n = plane.n
    end
    return isect
end

function orthobasis(n::Vect)
    b1 = zero(n)
    b2 = n
    if -0.6 < n.x < 0.6
        @set! b1.x = 1
    elseif -0.6 < n.y < 0.6
        @set! b1.y = 1;
    elseif -0.6 < n.z < 0.6
        @set! b1.z = 1;
    else
        @set! b1.x = 1
    end
    b0 = normalize(b1 × b2)
    b1 = normalize(b2 × b0)
    return b0, b1, b2
end

@fastmath @muladd function ambient_occlusion(isect::Intersect{T,V}, scene::State) where {T,V}
    @unpack nao_samples = scene
    ntheta = nao_samples
    nphi   = nao_samples
    ε = 0.0001
    p = isect.p + ε * isect.n
    b0, b1, b2 = orthobasis(isect.n)
    occlusion = zero(T)
    for j in 1:ntheta, i in 1:nphi
        theta = sqrt(rand())
        phi = 2pi * rand()

        x = cos(phi) * theta
        y = sin(phi) * theta
        z = sqrt(1 - theta^2)

        dir = x * b0 + y * b1 + z * b2
        ray = Ray(p, dir)

        occisect = Intersect(1e17, zero(V), zero(V), 0)
        occisect = ray_intersect(occisect, ray, scene.spheres[1])
        occisect = ray_intersect(occisect, ray, scene.spheres[2])
        occisect = ray_intersect(occisect, ray, scene.spheres[3])
        occisect = ray_intersect(occisect, ray, scene.plane)
        iszero(occisect.hit) || (occlusion += 1.0)
    end
    occlusion = (ntheta * nphi - occlusion) / (ntheta * nphi)
    return occlusion
end

function render!(img, h::Int, w::Int, nsubsamples::Int, scene::State)
    ss = inv(nsubsamples^2)
    fill!(img, 0)

    for y in 0:h-1, x in 0:w-1
        for v in 0:nsubsamples-1, u in 0:nsubsamples-1
            px = (x + (u / nsubsamples) - (w / 2.0)) / (w / 2.0)
            py = -(y + (v / nsubsamples) - (h / 2.0)) / (h / 2.0)
            dir = normalize(Vect(px, py, -1))
            ray = Ray(zero(dir), dir)
            isect = Intersect(1e17, zero(dir), zero(dir), 0)
            isect = ray_intersect(isect, ray, scene.spheres[1])
            isect = ray_intersect(isect, ray, scene.spheres[2])
            isect = ray_intersect(isect, ray, scene.spheres[3])
            isect = ray_intersect(isect, ray, scene.plane)

            if !iszero(isect.hit)
                col = ambient_occlusion(isect, scene)
                img[y+1, x+1] += Gray(col)
            end
        end

        img[y+1, x+1] *= ss
    end
    return img
end

function init_scene(;width=256, height=256, nsubsamples=2, nao_samples=8)
    spheres = (
               Sphere(Vect(-2.0, 0.0, -3.5), 0.5),
               Sphere(Vect(-0.5, 0.0, -3.0), 0.5),
               Sphere(Vect(1.0, 0.0, -2.2), 0.5),
              )
    plane = Plane(Vect(0.0, -0.5, 0.0), Vect(0.0, 1.0, 0.0))
    return State(spheres, plane, width, height, nsubsamples, nao_samples)
end

function main(;width=256, height=256, nsubsamples=2, nao_samples=8)
    scene = init_scene(;
                        width=width,
                        height=height,
                        nsubsamples=nsubsamples,
                        nao_samples=nao_samples,
                      )
    img = Array{RGB{Float64}}(undef, scene.height, scene.width)
    render!(img, scene.height, scene.width, scene.nsubsamples, scene)
    save(File(format"PNG", "ao.png"), img)
    return img
end

#=
function saveppm(fname::String, w::Int, h::Int, img::Array{UInt8,3})
    open(fname, "w+") do f
        println(f, "P6")
        println(f, "$w $h")
        println(f, "255")
        write(f, vec(UInt8.(img)))
    end
    return fname
end
=#

end # module
