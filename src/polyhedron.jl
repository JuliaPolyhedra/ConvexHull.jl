export ConvexHullLibrary

mutable struct ConvexHullLibrary <: PolyhedraLibrary
    precision::Symbol

    function ConvexHullLibrary(precision::Symbol=:float)
        if !(precision in [:float, :exact])
            error("Invalid precision, it should be :float or :exact")
        end
        new(precision)
    end
end

Polyhedra.similar_library(l::ConvexHullLibrary, ::FullDim, ::Type{<:AbstractFloat}) = ConvexHullLibrary(:float)
Polyhedra.similar_library(l::ConvexHullLibrary, ::FullDim, ::Type) = ConvexHullLibrary(:exact)

mutable struct ConvexHullPolyhedron{N, T} <: Polyhedron{N, T}
    ine::Nullable{HRepresentation{N, T}}
    inel::Nullable{LiftedHRepresentation{N, T}}
    ext::Nullable{VRepresentation{N, T}}
    extl::Nullable{LiftedVRepresentation{N, T}}
    noredundantinequality::Bool
    noredundantgenerator::Bool

    function ConvexHullPolyhedron{N, T}(ine, ext, nri::Bool, nrg::Bool) where {N, T}
        new(ine, Nullable{LiftedHRepresentation{N, T}}(), ext, Nullable{LiftedVRepresentation{N, T}}(), nri, nrg)
    end
    function ConvexHullPolyhedron{N, T}(ine::HRepresentation{N, T}) where {N, T}
        new(ine, Nullable{LiftedHRepresentation{N, T}}(), Nullable{VRepresentation{N, T}}(), Nullable{LiftedVRepresentation{N, T}}(), false, false)
    end
    function ConvexHullPolyhedron{N, T}(ext::VRepresentation{N, T}) where {N, T}
        new(Nullable{HRepresentation{N, T}}(), Nullable{LiftedHRepresentation{N, T}}(), ext, Nullable{LiftedVRepresentation{N, T}}(), false, false)
    end
end

Polyhedra.library(p::ConvexHullPolyhedron{N, <:AbstractFloat}) where N = ConvexHullLibrary(:float)
Polyhedra.library(p::ConvexHullPolyhedron) = ConvexHullLibrary(:exact)
Polyhedra.similar_type(::Type{<:ConvexHullPolyhedron}, ::FullDim{N}, ::Type{T}) where {N, T} = ConvexHullPolyhedron{N, T}

function Polyhedra.arraytype(p::ConvexHullPolyhedron)
    if isnull(p.ine) && !isnull(p.inel)
        p.ine = get(p.inel)
    end
    if isnull(p.ext) && !isnull(p.extl)
        p.ext = get(p.extl)
    end
    if isnull(p.ine)
        Polyhedra.arraytype(get(p.ext))
    elseif isnull(p.ext)
        Polyhedra.arraytype(get(p.ine))
    else
        @assert Polyhedra.arraytype(get(p.ine)) == Polyhedra.arraytype(get(p.ext))
        Polyhedra.arraytype(get(p.ine))
    end
end

# Helpers
function getine(p::ConvexHullPolyhedron)
    if isnull(p.ine)
        if !isnull(p.inel)
            p.ine = p.inel
        else
            p.ine = double_description(getextl(p))
            p.inel = nothing
            p.noredundantinequality = true
        end
    end
    get(p.ine)
end
function getinel(p::ConvexHullPolyhedron)
    if isnull(p.inel)
        p.inel = LiftedHRepresentation(getine(p))
    end
    get(p.inel)
end
function getext(p::ConvexHullPolyhedron)
    if isnull(p.ext)
        if !isnull(p.extl)
            p.ext = p.extl
        else
            p.ext = double_description(getinel(p))
            p.extl = nothing
            p.noredundantgenerator = true
        end
    end
    get(p.ext)
end
function getextl(p::ConvexHullPolyhedron)
    if isnull(p.extl)
        p.extl = LiftedVRepresentation(getext(p))
    end
    get(p.extl)
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
    if !isnull(p.ine)
        ine = copy(get(p.ine))
    end
    ext = nothing
    if !isnull(p.ext)
        ext = copy(get(p.ext))
    end
    ConvexHullPolyhedron{N, T}(ine, ext, p.noredundantinequality, p.noredundantgenerator)
end
function Polyhedra.hrepiscomputed(p::ConvexHullPolyhedron)
    !isnull(p.ine)
end
function Polyhedra.hrep(p::ConvexHullPolyhedron)
    getine(p)
end
function Polyhedra.vrepiscomputed(p::ConvexHullPolyhedron)
    !isnull(p.ext)
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
