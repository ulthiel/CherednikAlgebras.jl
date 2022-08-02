################################################################################
# Gaudin operators
################################################################################

export GaudinOperator

struct GaudinOperator
    rep::Representation{<:Gapjm.ComplexReflectionGroup, <:Oscar.Field}
    par::Parameter
    y::Vector{<:Oscar.FieldElem}

    matrix::Oscar.MatElem
    irrep_num::Int
end

function GaudinOperator(rep::Representation{<:Gapjm.ComplexReflectionGroup,<: Oscar.Field}, c::Parameter{ParameterTypes.EG}, y::Vector{<:Oscar.FieldElem})

    W = group(rep)
    n = dim(W)
    if length(y) != n
        throw(ArgumentError("Length of vector needs to be equal to dimension of group"))
    end
    refls = reflections(W)
    K = base_ring(c)
    𝕜 = codomain(c)
    S,y = PolynomialRing(𝕜, ["y"*string(i) for i=1:n])
    L = FractionField(S)
    y = map(L, y)
    α = coroots(W)
    α_K = [ matrix(K,1,n,α[i]) for i=1:length(α) if isassigned(α,i) ]

    return α_K

end
