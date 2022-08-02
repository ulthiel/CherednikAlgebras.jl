################################################################################
# Gapjm specific tools
################################################################################

export reflection_classes, Representation, position_reflection_class, hyperplane_orbits_pairs, hyperplane_orbits_pairs_full

################################################################################
# Conversions to OSCAR types
################################################################################

# Conversion of Gapjm.Cyc
function (K::QQFields)(x::Gapjm.Cyc)
    n = Gapjm.conductor(x)
    L,z = CyclotomicField(n)
    zK = K(z)
    xK = zero(K)
    for (i,c) in Gapjm.pairs(x)
        xK += c*zK^i
    end
    return xK
end

# Conversion of Gapjm.Root1
function (K::QQFields)(x::Gapjm.Root1)
    n = Gapjm.order(x)
    r = Gapjm.exponent(x)
    L,z = CyclotomicField(n)
    return K(z)^r
end

# Base field for (vectors of) matrices over ZZ/QQ/Cyc
function Gapjm_conductor(mats::Vector{Matrix{T}}) where T <: Union{Integer, Rational}
    return 1
end

function Gapjm_conductor(mats::Vector{Matrix{T}}) where T <: Gapjm.Cyc
    n = 1
    for A in mats
        nA = lcm([d for d in union([ Set([Gapjm.conductor(x) for x in Set(A)])])])
        n = lcm(n,nA)
    end
    return n
end

function Gapjm_conductor(A::Matrix{T}) where T <: Union{Integer, Rational, Gapjm.Cyc}
    return Gapjm_conductor([A])
end

function base_field(mats::Vector{Matrix{T}}) where T <: Union{Integer, Rational, Gapjm.Cyc}
    n = Gapjm_conductor(mats)
    if n == 1
        return QQ
    else
        K,z = CyclotomicField(n)
        return K
    end
end

function base_field(A::Matrix{T}) where T <: Union{Integer, Rational, Gapjm.Cyc}
    return base_field([A])
end

# Base field of complex reflection group
# reps=true takes all irreducible reps into account
function base_field(W::Gapjm.ComplexReflectionGroup ; reps=false)
    n = Gapjm_conductor(reflrep(W))
    if reps
        for rep in representations(W)
            n = lcm(n,Gapjm_conductor(rep))
        end
    end
    if n == 1
        return QQ
    else
        K,z = CyclotomicField(n)
        return K
    end
end

# Conversion of complex reflection group
function matrix_group(K::QQFields, W::Gapjm.ComplexReflectionGroup)
    gens = [matrix(K,A) for A in reflrep(W)]
    return matrix_group(gens)
end

function matrix_group(W::Gapjm.ComplexReflectionGroup)
    return matrix_group(base_field(W), W)
end

################################################################################
# Representations
################################################################################
struct Representation{T<:Gapjm.ComplexReflectionGroup,S<:Oscar.Field}
    group::T
    base_ring::S
    mats::Vector{<:Oscar.MatElem}
    dim::Int
end

function group(rho::Representation)
    return rho.group
end

function base_ring(rho::Representation)
    return rho.base_ring
end

function dim(rho::Representation)
    return rho.dim
end

function Base.show(io::IO, rho::Representation)
    print(io, "Representation of ", group(rho), " of dimension ", dim(rho), " over ", base_ring(rho) )
end

function Representation(W::Gapjm.ComplexReflectionGroup, K::QQFields, mats::Vector{Matrix{T}}) where T <: Union{Integer, Rational, Gapjm.Cyc}
    mats_K = [ matrix(K, m) for m in mats ]
    if any(m->is_square(m)==false, mats_K)
        throw(ArgumentError("Matrices must be square"))
    end
    d = ncols(mats_K[1])
    if any(m->ncols(m)!=d, mats)
        throw(ArgumentError("Matrices must be of the same size"))
    end
    return Representation{typeof(W),typeof(K)}(W, K, mats_K, d)
end

function Representation(W::Gapjm.ComplexReflectionGroup, mats::Vector{Matrix{T}}) where T <: Union{Integer, Rational, Gapjm.Cyc}
    K = base_field(mats)
    return Representation(W, K, mats)
end

function Representation(W::Gapjm.ComplexReflectionGroup, i::Int)
    return Representation(W, representation(W,i))
end

function (rho::Representation{T,S})(x::Vector{Int}) where T <: Gapjm.ComplexReflectionGroup where S <: Oscar.Field
    res = identity_matrix(base_field(rho), dim(rho))
    for i=1:length(x)
        res *= rho.mats[x[i]]
    end
    return res
end

function (rho::Representation{T,S})(i::Int) where T <: Gapjm.ComplexReflectionGroup where S <: Oscar.Field
    return rho.mats[i]
end

function (rho::Representation{T,S})(w::Gapjm.Perm) where T<:Gapjm.ComplexReflectionGroup where S<:Oscar.Field
    W = group(rho)
    w_word = Gapjm.word(W,w)
    return rho(w_word)
end

function (rho::Representation{T,S})(r::Gapjm.Reflection) where T<:Gapjm.ComplexReflectionGroup where S<:Oscar.Field
    return rho(word(r))
end


################################################################################
# Reflections et al
################################################################################
"""
    function reflection_classes(W::Gapjm.ComplexReflectionGroup)

The positions of the conjugacy classes of reflections in the list of classes of W. The list is ordered as in ``hyperplane_orbits``.
"""
function reflection_classes(W::Gapjm.ComplexReflectionGroup)
    if !haskey(W.prop, :reflection_classes)
        reflcl = Int[]
        𝒜 = hyperplane_orbits(W)
        for Ω ∈ 𝒜
            append!(reflcl, Ω.cl_s)
        end
        W.prop[:reflection_classes] = reflcl
    end
    return W.prop[:reflection_classes]
end

"""
    function position_reflection_class(W::Gapjm.ComplexReflectionGroup, w::Gapjm.Perm)

If w is a reflection, return the position of the conjugacy class of w in the
list of conjugacy classes of reflections.
"""
function position_reflection_class(W::Gapjm.ComplexReflectionGroup, w::Gapjm.Perm)
    return findfirst(isequal(position_class(W,w)), reflection_classes(W))
end

"""
    function hyperplane_orbits_pairs(W::Gapjm.ComplexReflectionGroup)

Let 𝒜 be the set of reflecting hyperplanes of W. This function returns the list of pairs (i,j) where 1 ≤ i ≤ #𝒜/W and 1 ≤ j ≤ e_Ω-1, where e_Ω is the order of the stabilizer W_H for any H ∈ Ω.
"""
function hyperplane_orbits_pairs(W::Gapjm.ComplexReflectionGroup)
    return [ (i,j) for i=1:length(hyperplane_orbits(W)) for j=1:hyperplane_orbits(W)[i].order-1 ]
end

"""
    function hyperplane_orbits_pairs_full(W::Gapjm.ComplexReflectionGroup)

Let 𝒜 be the set of reflecting hyperplanes of W. This function returns the list of pairs (i,j) where 1 ≤ i ≤ #𝒜/W and 0 ≤ j ≤ e_Ω-1, where e_Ω is the order of the stabilizer W_H for any H ∈ Ω.
"""
function hyperplane_orbits_pairs_full(W::Gapjm.ComplexReflectionGroup)
    return [ (i,j) for i=1:length(hyperplane_orbits(W)) for j=0:hyperplane_orbits(W)[i].order-1 ]
end

"""
    function dim(W::Gapjm.ComplexReflectionGroup)

Dimension of the reflection representation of W.
"""
function dim(W::Gapjm.ComplexReflectionGroup)
    return size(reflrep(W)[1])[1]
end
