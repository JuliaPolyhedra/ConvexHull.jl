export ConvexHullLib

type ConvexHullLib <: PolyhedraLibrary
    precision::Symbol

    function ConvexHullLib(precision::Symbol=:float)
        if !(precision in [:float, :exact])
            error("Invalid precision, it should be :float or :exact")
        end
        new(precision)
    end
end

type ConvexHullPolyhedron{N, T} <: Polyhedron{N, T}
    ine::Nullable{HRepresentation{N, T}}
    inel::Nullable{LiftedHRepresentation{N, T}}
    ext::Nullable{VRepresentation{N, T}}
    extl::Nullable{LiftedVRepresentation{N, T}}
    hlinearitydetected::Bool
    vlinearitydetected::Bool
    noredundantinequality::Bool
    noredundantgenerator::Bool

    function ConvexHullPolyhedron(ine::HRepresentation{N, T}, ext::VRepresentation{N, T}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool)
        new(ine, Nullable{LiftedHRepresentation{N, T}}(), ext, Nullable{LiftedVRepresentation{N, T}}(), hld, vld, nri, nrg)
    end
    function ConvexHullPolyhedron(ine::HRepresentation{N, T})
        new(ine, Nullable{LiftedHRepresentation{N, T}}(), Nullable{VRepresentation{N, T}}(), Nullable{LiftedVRepresentation{N, T}}(), false, false, false, false)
    end
    function ConvexHullPolyhedron(ext::VRepresentation{N, T})
        new(Nullable{HRepresentation{N, T}}(), Nullable{LiftedHRepresentation{N, T}}(), ext, Nullable{LiftedVRepresentation{N, T}}(), false, false, false, false)
    end
end

# ine may decompose fast but if ine is nothing I do not want to ask to compute it to see the type it is
# saying false normally do not give troubles
decomposedhfast{N, T}(::Type{ConvexHullPolyhedron{N, T}}) = false
decomposedvfast{N, T}(::Type{ConvexHullPolyhedron{N, T}}) = false
decomposedhfast{N, T}(::ConvexHullPolyhedron{N, T}) = decomposedhfast(ConvexHullPolyhedron{N, T})
decomposedvfast{N, T}(::ConvexHullPolyhedron{N, T}) = decomposedvfast(ConvexHullPolyhedron{N, T})

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
    hlinearitydetected = false
    vlinearitydetected = false
    noredundantinequality = false
    noredundantgenerator = false
end
function updateine!{N,T}(p::ConvexHullPolyhedron{N}, ine::HRepresentation{N, T})
    clearfield!(p)
    p.ine = ine
end
function updateext!{N,T}(p::ConvexHullPolyhedron{N}, ext::VRepresentation{N, T})
    clearfield!(p)
    p.ext = ext
end


# Implementation of Polyhedron's mandatory interface
function polytypeforprecision(precision::Symbol)
  if !(precision in (:float, :exact))
    error("precision should be :float or :exact, you gave $precision")
  end
  precision == :float ? Float64 : Rational{BigInt}
end

function polyhedron{N}(repit::Union{Representation{N},HRepIterator{N},VRepIterator{N}}, lib::ConvexHullLib)
    T = polytypeforprecision(lib.precision)
    ConvexHullPolyhedron{N, T}(repit)
end

getlibraryfor{T<:Real}(::ConvexHullPolyhedron, n::Int, ::Type{T}) = ConvexHullLib(:exact)
getlibraryfor{T<:AbstractFloat}(::ConvexHullPolyhedron, n::Int, ::Type{T}) = ConvexHullLib(:float)

(::Type{ConvexHullPolyhedron{N, T}}){N, T}(it::HRepIterator{N,T}) = ConvexHullPolyhedron{N, T}(LiftedHRepresentation{N,T}(it))
(::Type{ConvexHullPolyhedron{N, T}}){N, T}(it::VRepIterator{N,T}) = ConvexHullPolyhedron{N, T}(LiftedVRepresentation{N,T}(it))

function (::Type{ConvexHullPolyhedron{N, T}}){N, T}(; eqs=nothing, ineqs=nothing, points=nothing, rays=nothing)
    noth = eqs === nothing && ineqs === nothing
    notv = points === nothing && rays === nothing
    if noth && notv
        error("ConvexHullPolyhedron should have at least one iterator to be built")
    end
    if !noth && !notv
        error("ConvexHullPolyhedron constructed with a combination of eqs/ineqs with points/rays")
    end
    if notv
        ConvexHullPolyhedron{N, T}(LiftedHRepresentation{N, T}(eqs=eqs, ineqs=ineqs))
    else
        ConvexHullPolyhedron{N, T}(LiftedVRepresentation{N, T}(points=points, rays=rays))
    end
end

function Base.copy{N}(p::ConvexHullPolyhedron{N})
    ine = nothing
    if !isnull(p.ine)
        ine = copy(get(p.ine))
    end
    ext = nothing
    if !isnull(p.ext)
        ext = copy(get(p.ext))
    end
    ConvexHullPolyhedron{N}(ine, ext, p.hlinearitydetected, p.vlinearitydetected, p.noredundantinequality, p.noredundantgenerator)
end
function Base.push!{N, T}(p::ConvexHullPolyhedron{N, T}, ine::HRepresentation{N})
    updateine!(p, intersect(getine(p), changeeltype(ine, T)))
end
function Base.push!{N, T}(p::ConvexHullPolyhedron{N, T}, ext::VRepresentation{N})
    updateext!(p, getext(p) + changeeltype(ext, T))
end
function hrepiscomputed(p::ConvexHullPolyhedron)
    !isnull(p.ine)
end
function hrep(p::ConvexHullPolyhedron)
    copy(getine(p))
end
function vrepiscomputed(p::ConvexHullPolyhedron)
    !isnull(p.ext)
end
function vrep(p::ConvexHullPolyhedron)
    copy(getext(p))
end
function detecthlinearities!(p::ConvexHullPolyhedron)
    warn("detecthlinearities! not supported yet")
end
function detectvlinearities!(p::ConvexHullPolyhedron)
    warn("detectvlinearities! not supported yet")
end
function removehredundancy!(p::ConvexHullPolyhedron)
    warn("removehredundancy! not supported yet")
end
function removevredundancy!(p::ConvexHullPolyhedron)
    warn("removevredundancy! not supported yet")
end

for f in [:hashreps, :nhreps, :starthrep, :hasineqs, :nineqs, :startineq, :haseqs, :neqs, :starteq]
    @eval $f(p::ConvexHullPolyhedron) = $f(getine(p))
end
for f in [:donehrep, :nexthrep, :doneineq, :nextineq, :doneeq, :nexteq]
    @eval $f(p::ConvexHullPolyhedron, state) = $f(getine(p), state)
end

for f in [:hasvreps, :nvreps, :startvrep, :haspoints, :npoints, :startpoint, :hasrays, :nrays, :startray]
    @eval $f(p::ConvexHullPolyhedron) = $f(getext(p))
end
for f in [:donevrep, :nextvrep, :donepoint, :nextpoint, :doneray, :nextray]
    @eval $f(p::ConvexHullPolyhedron, state) = $f(getext(p), state)
end
