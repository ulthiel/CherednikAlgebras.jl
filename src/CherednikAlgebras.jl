################################################################################
# A Julia package to work with rational Chererednik algebras et al.
#
# This package uses:
#   * OSCAR (https://github.com/oscar-system/Oscar.jl) for general computer
#     algebra.
#   * Gapjm (https://github.com/jmichel7/Gapjm.jl) for complex reflection groups
#     et al.
#
# By Ulrich Thiel (University of Kaiserslautern), 2022
################################################################################

module CherednikAlgebras

################################################################################
# Imports
################################################################################
import Oscar:
    Oscar, matrix_group, base_field, gens, ZZ, QQ, CyclotomicField, QQBar,
    abelian_closure, matrix, gen, gens, base_ring, dim, parent, domain,
    codomain, is_rational, hom, PolynomialRing, FractionField, quo, ngens, is_square, ncols, identity_matrix, kernel, is_zero, zero_matrix, ideal, nrows

import Gapjm:
    Gapjm, crg, hyperplane_orbits, reflrep, representations, position_class, nref, nhyp, CharTable, Diagram, charinfo, charnames, classnames, classinfo,
    representation, reflections, coroots, roots

################################################################################
# Define global QQAb
# Have to do this inside the __init__ function, following
# https://nemocas.github.io/Nemo.jl/dev/misc/#Global-variables-and-precompilation
################################################################################
function __init__()
  global QQAb, QQAbgen = abelian_closure(QQ)
end

################################################################################
# Exports from above (for convenience)
################################################################################

# From OSCAR
export
    Oscar, matrix_group, base_field, gens, ZZ, QQ, CyclotomicField, QQBar,
    abelian_closure, matrix, gen, gens, base_ring, dim, parent, domain,
    codomain, is_rational, hom, PolynomialRing, FractionField, quo, ngens, ideal

# From Gapjm
export
    crg, hyperplane_orbits, reflrep, representations, position_class, nref, nhyp, CharTable, Diagram, charinfo, charnames, classnames, classinfo, representation, reflections, coroots, roots

# Global QQAb
export QQAb, QQAbgen

################################################################################
# Union types (for convenience)
################################################################################

# Union type for OSCAR algebraic field extensions of Q
QQFields = Union{Oscar.FlintRationalField, Oscar.AnticNumberField, Oscar.QQAbField, Oscar.CalciumQQBarField}

# Union type for elements of QQFields
QQFieldsElem = Union{Int, Rational, Oscar.fmpz, Oscar.fmpq, Oscar.nf_elem, Oscar.qqbar, Oscar.QQAbElem}

################################################################################
# Includes
################################################################################
include("Gapjm_tools.jl")
include("Parameters.jl")
include("GaudinOperators.jl")

#
# include("SymplecticDoubling.jl")

end # module
