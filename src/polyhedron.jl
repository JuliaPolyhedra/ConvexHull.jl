export ConvexHullLibrary

mutable struct ConvexHullLibrary <: Polyhedra.Library
    precision::Symbol

    function ConvexHullLibrary(precision::Symbol=:float)
        if !(precision in [:float, :exact])
            error("Invalid precision, it should be :float or :exact")
        end
        new(precision)
    end
end

Polyhedra.similar_library(l::ConvexHullLibrary, ::Polyhedra.FullDim, ::Type{<:AbstractFloat}) = ConvexHullLibrary(:float)
Polyhedra.similar_library(l::ConvexHullLibrary, ::Polyhedra.FullDim, ::Type) = ConvexHullLibrary(:exact)

mutable struct ConvexHullPolyhedron{T} <: Polyhedron{T}
    ine::Union{HRepresentation{T}, Nothing}
    inel::Union{LiftedHRepresentation{T}, Nothing}
    ext::Union{VRepresentation{T}, Nothing}
    extl::Union{LiftedVRepresentation{T}, Nothing}
    noredundantinequality::Bool
    noredundantgenerator::Bool

    function ConvexHullPolyhedron{T}(ine, ext, nri::Bool, nrg::Bool) where {T}
        new(ine, nothing, ext, nothing, nri, nrg)
    end
    function ConvexHullPolyhedron{T}(ine::HRepresentation{T}) where {T}
        new(ine, nothing, nothing, nothing, false, false)
    end
    function ConvexHullPolyhedron{T}(ext::VRepresentation{T}) where {T}
        new(nothing, nothing, ext, nothing, false, false)
    end
end

Polyhedra.library(p::ConvexHullPolyhedron{<:AbstractFloat}) = ConvexHullLibrary(:float)
Polyhedra.library(p::ConvexHullPolyhedron) = ConvexHullLibrary(:exact)
Polyhedra.similar_type(::Type{<:ConvexHullPolyhedron}, ::Polyhedra.FullDim, ::Type{T}) where {T} = ConvexHullPolyhedron{T}

function Polyhedra.vectortype(p::ConvexHullPolyhedron)
    if isnothing(p.ine) && !isnothing(p.inel)
        p.ine = p.inel
    end
    if isnothing(p.ext) && !isnothing(p.extl)
        p.ext = p.extl
    end
    if isnothing(p.ine)
        Polyhedra.vectortype(p.ext)
    elseif isnothing(p.ext)
        Polyhedra.vectortype(p.ine)
    else
        @assert Polyhedra.vectortype(p.ine) == Polyhedra.vectortype(p.ext)
        Polyhedra.vectortype(p.ine)
    end
end

# Helpers
function getine(p::ConvexHullPolyhedron)
    if isnothing(p.ine)
        if !isnothing(p.inel)
            p.ine = p.inel
        else
            p.ine = double_description(getextl(p))
            p.inel = nothing
            p.noredundantinequality = true
        end
    end
    p.ine
end
function getinel(p::ConvexHullPolyhedron)
    if isnothing(p.inel)
        p.inel = LiftedHRepresentation(getine(p))
    end
    p.inel
end
function getext(p::ConvexHullPolyhedron)
    if isnothing(p.ext)
        if !isnothing(p.extl)
            p.ext = p.extl
        else
            p.ext = double_description(getinel(p))
            p.extl = nothing
            p.noredundantgenerator = true
        end
    end
    p.ext
end
function getextl(p::ConvexHullPolyhedron)
    if isnothing(p.extl)
        p.extl = LiftedVRepresentation(getext(p))
    end
    p.extl
end

function clearfield!(p::ConvexHullPolyhedron)
    p.ine = nothing
    p.inel = nothing
    p.ext = nothing
    p.extl = nothing
    p.noredundantinequality = false
    p.noredundantgenerator = false
end

# Implementation of Polyhedron's mandatory interface
function polytypeforprecision(precision::Symbol)
  if !(precision in (:float, :exact))
    error("precision should be :float or :exact, you gave $precision")
  end
  precision == :float ? Float64 : Rational{BigInt}
end

function Polyhedra.polyhedron(rep::Representation{N}, lib::ConvexHullLibrary) where N
    T = polytypeforprecision(lib.precision)
    ConvexHullPolyhedron{N, T}(rep)
end

ConvexHullPolyhedron{N, T}(it::Polyhedra.HIt{N}...) where {N, T} = ConvexHullPolyhedron{N, T}(LiftedHRepresentation{N, T}(it...))
ConvexHullPolyhedron{N, T}(it::Polyhedra.VIt{N}...) where {N, T} = ConvexHullPolyhedron{N, T}(LiftedVRepresentation{N, T}(it...))

function Base.copy(p::ConvexHullPolyhedron{N, T}) where {N, T}
    ine = nothing
    if !isnothing(p.ine)
        ine = copy(get(p.ine))
    end
    ext = nothing
    if !isnothing(p.ext)
        ext = copy(get(p.ext))
    end
    ConvexHullPolyhedron{N, T}(ine, ext, p.noredundantinequality, p.noredundantgenerator)
end
function Polyhedra.hrepiscomputed(p::ConvexHullPolyhedron)
    !isnothing(p.ine)
end
function Polyhedra.hrep(p::ConvexHullPolyhedron)
    getine(p)
end
function Polyhedra.vrepiscomputed(p::ConvexHullPolyhedron)
    !isnothing(p.ext)
end
function Polyhedra.vrep(p::ConvexHullPolyhedron)
    getext(p)
end

function Polyhedra.sethrep!(p::ConvexHullPolyhedron{N}, h::HRepresentation{N}) where N
    p.ine = h
    p.inel = nothing
end
function Polyhedra.setvrep!(p::ConvexHullPolyhedron{N}, v::VRepresentation{N}) where N
    p.ext = v
    p.extl = nothing
end
function resethrep!(p::ConvexHullPolyhedron{N}, h::HRepresentation{N}) where N
    clearfield!(p)
    p.ine = h
end
function resetvrep!(p::ConvexHullPolyhedron{N}, v::VRepresentation{N}) where N
    clearfield!(p)
    p.ext = v
end
