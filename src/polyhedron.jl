mutable struct Library <: Polyhedra.Library
    precision::Symbol
    solver::Polyhedra.SolverOrNot

    function Library(precision::Symbol=:float, solver=nothing)
        if !(precision in [:float, :exact])
            error("Invalid precision, it should be :float or :exact")
        end
        new(precision, solver)
    end
end

Polyhedra.similar_library(l::Library, ::Polyhedra.FullDim, ::Type{<:AbstractFloat}) = Library(:float)
Polyhedra.similar_library(l::Library, ::Polyhedra.FullDim, ::Type) = Library(:exact)

const HRepT{T} = LiftedHRepresentation{T, Matrix{T}}
const VRepT{T} = LiftedVRepresentation{T, Matrix{T}}

mutable struct Polyhedron{T} <: Polyhedra.Polyhedron{T}
    hrep::Union{HRepT{T}, Nothing}
    vrep::Union{VRepT{T}, Nothing}
    noredundantinequality::Bool
    noredundantgenerator::Bool
    solver::Polyhedra.SolverOrNot
end

function Polyhedron{T}(hrep::HRepresentation, solver::Polyhedra.SolverOrNot) where {T}
    Polyhedron{T}(convert(HRepT{T}, hrep), nothing, false, false, solver)
end
function Polyhedron{T}(vrep::VRepresentation, solver::Polyhedra.SolverOrNot) where {T}
    Polyhedron{T}(nothing, convert(VRepT{T}, vrep), false, false, solver)
end
#function Polyhedron{T}(rep::Representation, solver::Polyhedra.SolverOrNot) where {T}
#    return Polyhedron{T}(Polyhedra.change_coefficient_type(rep, T), solver)
#end

Polyhedra.FullDim(p::Polyhedron) = Polyhedra.FullDim_rep(p.hrep, p.vrep)
Polyhedra.library(p::Polyhedron{<:AbstractFloat}) = Library(:float, p.solver)
Polyhedra.library(p::Polyhedron) = Library(:exact, p.solver)
Polyhedra.default_solver(p::Polyhedron; T = nothing) = p.solver
Polyhedra.supportssolver(::Type{<:Polyhedron}) = true

function Polyhedra.hvectortype(::Type{Polyhedron{T}}) where {T}
    return Polyhedra.hvectortype(HRepT{T})
end
function Polyhedra.vvectortype(::Type{Polyhedron{T}}) where {T}
    return Polyhedra.vvectortype(VRepT{T})
end

Polyhedra.similar_type(::Type{<:Polyhedron}, ::Polyhedra.FullDim, ::Type{T}) where {T} = Polyhedron{T}

function Polyhedron{T}(d::Polyhedra.FullDim, it::Polyhedra.HIt...; solver=nothing) where {T}
    return Polyhedron{T}(HRepT{T}(d, it...), solver)
end
function Polyhedron{T}(d::Polyhedra.FullDim, it::Polyhedra.VIt...; solver=nothing) where {T}
    return Polyhedron{T}(VRepT{T}(d, it...), solver)
end

function clearfield!(p::Polyhedron)
    p.hrep = nothing
    p.vrep = nothing
    p.noredundantinequality = false
    p.noredundantgenerator = false
end

# Implementation of Polyhedron's mandatory interface
function polytypeforprecision(precision::Symbol)
  if !(precision in (:float, :exact))
    error("precision should be :float or :exact, you gave $precision")
  end
  return precision == :float ? Float64 : Rational{BigInt}
end

function Polyhedra.polyhedron(rep::Representation, lib::Library)
    T = polytypeforprecision(lib.precision)
    return Polyhedron{T}(rep, lib.solver)
end

function Base.copy(p::Polyhedron{T}) where {T}
    hrep = nothing
    if !isnothing(p.hrep)
        hrep = copy(p.hrep)
    end
    vrep = nothing
    if !isnothing(p.vrep)
        vrep = copy(p.vrep)
    end
    return Polyhedron{T}(hrep, vrep, p.noredundantinequality,
                         p.noredundantgenerator, p.solver)
end

function Polyhedra.hrepiscomputed(p::Polyhedron)
    !isnothing(p.hrep)
end
function Polyhedra.computehrep!(p::Polyhedron)
    # vrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.hrep = double_description(p.vrep)
    p.noredundantinequality = true
    return
end
function Polyhedra.hrep(p::Polyhedron)
    if !Polyhedra.hrepiscomputed(p)
        Polyhedra.computehrep!(p)
    end
    return p.hrep
end
function Polyhedra.vrepiscomputed(p::Polyhedron)
    !isnothing(p.vrep)
end
function Polyhedra.computevrep!(p::Polyhedron)
    # hrep(p) could trigger an infinite loop if both vrep and hrep are null
    p.vrep = double_description(p.hrep)
    p.noredundantgenerator = true
    return
end
function Polyhedra.vrep(p::Polyhedron)
    if !Polyhedra.vrepiscomputed(p)
        Polyhedra.computevrep!(p)
    end
    return p.vrep
end

function Polyhedra.sethrep!(p::Polyhedron{T}, h::HRepresentation) where T
    p.hrep = convert(HRepT{T}, h)
end
function Polyhedra.setvrep!(p::Polyhedron{T}, v::VRepresentation) where T
    p.vrep = convert(VRepT{T}, v)
end
function Polyhedra.resethrep!(p::Polyhedron, h::HRepresentation)
    clearfield!(p)
    Polyhedra.sethrep!(p, h)
end
function Polyhedra.resetvrep!(p::Polyhedron, v::VRepresentation)
    clearfield!(p)
    Polyhedra.setvrep!(p, v)
end
