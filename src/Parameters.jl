################################################################################
# Parameters for rational Cherednik algebras ################################################################################

################################################################################
# We create types for all the different types of parameters.
# These come into the submodule ParameterTypes so that they have their own
# namespace.
################################################################################
module ParameterTypes

# Imports from parent module
import ..Gapjm, ..reflection_classes, ..hyperplane_orbits

abstract type Abs end # Abstract supertype

# Generic functions
function group(par::Abs)
    return par.W
end

function type_str(par::T) where T <: Abs
    return string(nameof(T))
end

function Base.show(io::IO, par::T) where T <: Abs
    print(io, "Parameter of type ", type_str(par), " for ", group(par))
end

# Generic function to create parameter type for a group using a type string
function create(W::Gapjm.ComplexReflectionGroup, par_type_str::String)
    f = Symbol(par_type_str)
    return getfield(ParameterTypes, f)(W)
end

# Now, the various types

# Etingof-Ginzburg
struct EG <: Abs
    W::Gapjm.ComplexReflectionGroup
end

# Bellamy-Schedler-Thiel (κ)
struct BST <: Abs
    W::Gapjm.ComplexReflectionGroup
end

# McKay (k from BST)
struct McKay <: Abs
    W::Gapjm.ComplexReflectionGroup
end

# Ginzburg-Guay-Opdam-Rouquier
struct GGOR <: Abs
    W::Gapjm.ComplexReflectionGroup
end

# Bonnafé-Rouquier (C)
struct BR_C <: Abs
    W::Gapjm.ComplexReflectionGroup
end

# Bonnafé-Rouquier (K)
struct BR_K <: Abs
    W::Gapjm.ComplexReflectionGroup
end

# Broué-Malle-Rouquier
struct BMR <: Abs
    W::Gapjm.ComplexReflectionGroup
end

end # of module ParameterTypes

################################################################################
# Parameter spaces
################################################################################
export ParameterSpace, group, coordinate_ring, par_type

mutable struct ParameterSpace{T<:ParameterTypes.Abs}
    par_type::T
    coordinate_ring::Oscar.MPolyRing
    transformation_matrix::Dict{ParameterSpace, Oscar.MatElem}
end

# Generic functions
function par_type(P::ParameterSpace)
    return P.par_type
end

function coordinate_ring(P::ParameterSpace)
    return P.coordinate_ring
end

function group(P::ParameterSpace)
    return ParameterTypes.group(par_type(P))
end

function base_ring(P::ParameterSpace)
    return base_ring(coordinate_ring(P))
end

function dim(P::ParameterSpace)
    return ngens(coordinate_ring(P))
end

# Generic function to get transformation_matrix if inverse is available
function transformation_matrix(P::ParameterSpace, Q::ParameterSpace)

    if haskey(P.transformation_matrix, Q)
        return P.transformation_matrix[Q]
    end

    if hasmethod(transformation_matrix, (typeof(Q), typeof(P)))
        A = transformation_matrix(Q, P)
        P.transformation_matrix[Q] = inv(A)
        return P.transformation_matrix[Q]
    end

    throw(ArgumentError("Transformation matrix not implemented for this type"))

end

# Printing
function Base.show(io::IO, P::ParameterSpace)
    print(io, "Parameter space of type ", ParameterTypes.type_str(par_type(P)), " for ", group(P), " over ", base_ring(P) )
end

# Constructors
function ParameterSpace{T}(par_type::ParameterTypes.Abs, R::Oscar.MPolyRing) where T <: ParameterTypes.Abs
    return ParameterSpace{T}(par_type, R, Dict{ParameterSpace, Oscar.Matrix}())
end

function ParameterSpace(W::Gapjm.ComplexReflectionGroup, par_type_str::String, K::QQFields; kwargs...)
    par_type = ParameterTypes.create(W, par_type_str)
    return ParameterSpace(par_type, K; kwargs...)
end

function ParameterSpace(W::Gapjm.ComplexReflectionGroup, K::QQFields)
    return ParameterSpace(par_type, "EG", K)
end

function ParameterSpace(W::Gapjm.ComplexReflectionGroup, par_type_str::String; kwargs...)
    return ParameterSpace(W, par_type_str, QQAb; kwargs...)
end

function ParameterSpace(W::Gapjm.ComplexReflectionGroup; kwargs...)
    return ParameterSpace(W, "EG"; kwargs...)
end

function ParameterSpace(par_type::ParameterTypes.Abs)
    return ParameterSpace(par_type, QQAb)
end

# EG
function ParameterSpace(par_type::ParameterTypes.EG, K::QQFields; label="c")

    W = ParameterTypes.group(par_type)
    vars_labels = [ label*"("*string(i)*")" for i in 1:length(reflection_classes(W)) ]
    R, vars = PolynomialRing(K, vars_labels)
    return ParameterSpace{ParameterTypes.EG}(par_type, R)

end

# BST
function ParameterSpace(par_type::ParameterTypes.BST, K::QQFields; label="κ")

    W = ParameterTypes.group(par_type)
    vars_labels = [ label*"("*string(x[1])*","*string(x[2])*")" for x in hyperplane_orbits_pairs(W) ]
    R, vars = PolynomialRing(K, vars_labels)
    return ParameterSpace{ParameterTypes.BST}(par_type, R)

end

# GGOR
function ParameterSpace(par_type::ParameterTypes.GGOR, K::QQFields; label="k")

    W = ParameterTypes.group(par_type)
    vars_labels = [ label*"("*string(x[1])*","*string(x[2])*")" for x in hyperplane_orbits_pairs(W) ]
    R, vars = PolynomialRing(K, vars_labels)
    return ParameterSpace{ParameterTypes.GGOR}(par_type, R)

end

# BMR
function ParameterSpace(par_type::ParameterTypes.BMR, K::QQFields; label="u")

    W = ParameterTypes.group(par_type)
    vars_labels = [ label*"("*string(x[1])*"_"*string(x[2])*")" for x in hyperplane_orbits_pairs_full(W) ]
    R, vars = PolynomialRing(K, vars_labels)
    return ParameterSpace{ParameterTypes.BMR}(par_type, R)

end

################################################################################
# Parameters
################################################################################
export Parameter, generic_point, make_rational, is_generic, space, alg_hom, is_rational, parent, is_closed

struct Parameter{T<:ParameterTypes.Abs}
    parent::ParameterSpace{T}
    alg_hom::Oscar.MPolyAnyMap
end

# Generic functions
function parent(par::Parameter)
    return par.parent
end

function alg_hom(par::Parameter)
    return par.alg_hom
end

function domain(par::Parameter)
    return coordinate_ring(parent(par))
end

function codomain(par::Parameter)
    return codomain(alg_hom(par))
end

function par_type(par::Parameter)
    return par_type(parent(par))
end

function group(par::Parameter)
    return group(parent(par))
end

function base_ring(par::Parameter)
    return base_ring(parent(par))
end

function Base.show(io::IO, par::Parameter)
    P = parent(par)
    print(io, "Parameter of type ", ParameterTypes.type_str(par_type(P)), " for ", group(P), " with values ", "["*join([alg_hom(par)(x) for x in gens(domain(par))],", ")*"]")
end

function is_generic(par::Parameter{T}) where T <: ParameterTypes.Abs
    return is_zero(kernel(alg_hom(par)))
end

function is_rational(par::Parameter{T}) where T <: ParameterTypes.Abs
    𝕜 = codomain(par)
    return isa(𝕜, Oscar.Field)
end

function is_closed(par::Parameter{T}) where T <: ParameterTypes.Abs
    P = parent(par)
    R = coordinate_ring(P)
    K = base_ring(P)
    φ = alg_hom(par)

    # If φ maps into K then the kernel is maximal.
    # Otherwise we have to test.
    if codomain(φ) == K
        return true
    else
        I = kernel(φ)
        return dim(I) == 0
    end
end

function make_rational(par::Parameter{T}) where T <: ParameterTypes.Abs
    if is_rational(par)
        return par
    end
    𝕜 = codomain(par)
    𝕂 = FractionField(𝕜)
    P = parent(par)
    R = domain(par)
    φ = alg_hom(par)
    rat_φ = hom(R,𝕂,[𝕂(φ(x)) for x in gens(R)])
    return Parameter{T}(P, rat_φ)
end

# Constructors
function generic_point(P::ParameterSpace{T}) where T <: ParameterTypes.Abs
    R = coordinate_ring(P)
    φ = hom(R, R, gens(R))
    return Parameter{T}(P, φ)
end

function (P::ParameterSpace{T})() where T <: ParameterTypes.Abs
    return generic_point(P)
end

function (P::ParameterSpace{T})(pt::Oscar.MatElem{S}) where T <: ParameterTypes.Abs where S <: Oscar.RingElem

    if nrows(pt) != 1
        throw(ArgumentError("Matrix needs to have one row only"))
    end

    if ncols(pt) != dim(P)
        throw(ArgumentError("Matrix needs to have " * string(dim(P)) * " columns"))
    end

    R = coordinate_ring(P)
    𝕜 = base_ring(pt)
    φ = hom(R, 𝕜, vec(Array(pt)))

    return Parameter{T}(P, φ)
end

function (P::ParameterSpace{T})(pt::Vector{S}) where T <: ParameterTypes.Abs where S <: QQFieldsElem
    if length(pt) != ngens(coordinate_ring(P))
        throw(ArgumentError("Point needs to have " * string(dim(P)) * " coordinates"))
    end

    # We will coerce all values into the base ring of the space
    𝕜 = base_ring(P)
    vals = [𝕜(x) for x in pt]

    R = coordinate_ring(P)
    φ = hom(R, 𝕜, vals)

    return Parameter{T}(P, φ)
end

function (P::ParameterSpace{T})(I::Oscar.MPolyIdeal) where T <: ParameterTypes.Abs

    R = coordinate_ring(P)

    if base_ring(I) != R
        throw(ArgumentError("Ideal does not belong to the coordinate ring of the space"))
    end

    𝕜,q1 = quo(R,I)
    q = hom(R, 𝕜, [q1(x) for x in gens(R)])

    return Parameter{T}(P, q)
end

# Using transformation matrix
function (P::ParameterSpace)(par::Parameter)
    Q = parent(par)
    A = transformation_matrix(Q,P)
    φ = alg_hom(par)
    R = domain(φ)
    𝕜 = codomain(φ)
    v = matrix(𝕜, 1, ngens(R), [φ(x) for x in gens(R)])
    return P(v*A)
end

################################################################################
# Special functions for EG parameters ################################################################################
function (c::Parameter{ParameterTypes.EG})(i::Int)
    φ = alg_hom(c)
    return φ(gens(domain(φ))[i])
end

function (c::Parameter{ParameterTypes.EG})(w::Gapjm.Perm)

    # Get the position in cl of the cojugacy class number of w
    W = group(c)
    i = findfirst(isequal(position_class(W,w)), par_type(c).reflection_classes)
    if isnothing(i)
        throw(ArgumentError("Element is not a reflection"))
    end
    return c(i)

end

function (c::Parameter{ParameterTypes.EG})(r::Gapjm.Reflection)
    i = findfirst(isequal(position_class(r)), par_type(c).reflection_classes)
    return c(i)
end

################################################################################
# Special functions for BST parameters ################################################################################
function (κ::Parameter{ParameterTypes.BST})(i::Int, j::Int)
    W = group(κ)
    𝒜 = hyperplane_orbits(W)
    𝒜_labels = hyperplane_orbits_pairs(W)
    Ω = 𝒜[i]
    e_Ω = Ω.order
    φ = alg_hom(κ)

    #Take index within orbit Ω mod e_Ω
    j = mod(j,e_Ω)

    if j > 0
        l = findfirst(isequal((i,j)), 𝒜_labels)
        if isnothing(l)
            throw(ArgumentError("Index out of range"))
        end
        return φ(gens(domain(φ))[l])
    else
        # For j=0 return -Σ of the k_{i,l} in the orbit
        # (this is a relation that we have for BST parameters)
        v = zero(codomain(φ))
        for l=1:e_Ω-1
            v = v - κ(i,l)
        end
        return v
    end
end

################################################################################
# Special functions for GGOR parameters ################################################################################
function (k::Parameter{ParameterTypes.GGOR})(i::Int, j::Int)
    W = group(k)
    𝒜 = hyperplane_orbits(W)
    𝒜_labels = hyperplane_orbits_pairs(W)
    Ω = 𝒜[i]
    e_Ω = Ω.order
    φ = alg_hom(k)

    #Take index within orbit Ω mod e_Ω
    j = mod(j,e_Ω)

    if j > 0
        l = findfirst(isequal((i,j)), 𝒜_labels)
        if isnothing(l)
            throw(ArgumentError("Index out of range"))
        end
        return φ(gens(domain(φ))[l])
    else
        # For j=0 return 0
        # (this is a relation that we have for GGOR parameters)
        v = zero(codomain(φ))
        return v
    end
end

################################################################################
# Special functions for BMR parameters ################################################################################
function (u::Parameter{ParameterTypes.BMR})(i::Int, j::Int)
    W = group(u)
    𝒜 = hyperplane_orbits(W)
    𝒜_labels = hyperplane_orbits_pairs_full(W)
    Ω = 𝒜[i]
    e_Ω = Ω.order
    φ = alg_hom(u)

    #Take index within orbit Ω mod e_Ω
    j = mod(j,e_Ω)

    l = findfirst(isequal((i,j)), 𝒜_labels)
    if isnothing(l)
        throw(ArgumentError("Index out of range"))
    end
    return φ(gens(domain(φ))[l])
end

################################################################################
# Converting between from BST to EG
################################################################################
export transformation_matrix

function transformation_matrix(𝒦::ParameterSpace{ParameterTypes.BST}, 𝒞::ParameterSpace{ParameterTypes.EG})

    if haskey(𝒦.transformation_matrix, 𝒞)
        return 𝒦.transformation_matrix[𝒞]
    end

    if group(𝒞) != group(𝒦)
        throw(ArgumentError("Parameter spaces do not have the same group"))
    end

    if base_ring(𝒞) != base_ring(𝒦)
        throw(ArgumentError("Parameter spaces do not have the same base ring"))
    end

    W = group(𝒞)
    K = base_ring(𝒞)
    d = dim(𝒞)
    A = zero_matrix(K, d, d)
    l = 1 #row counter
    for Ω in hyperplane_orbits(W)
        e_Ω = Ω.order
        det = K(QQAbgen(e_Ω))

        #Translate class indices to reflection class indices
        refl_cl = [findfirst(isequal(i), reflection_classes(W)) for i in Ω.cl_s]

        for i=1:e_Ω-1
            for j in refl_cl
                A[l,j] = det^((i-1)*j) - det^(i*j)
            end
            l = l + 1
        end
    end

    𝒦.transformation_matrix[𝒞] = A

    return A
end
