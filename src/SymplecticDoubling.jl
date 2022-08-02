################################################################################
# Symplectic doubling of a matrix group.
################################################################################

export
    SymplecticDoubling,
    symplectic_doubling


struct SymplecticDoubling
    org_group::MatrixGroup
    group::MatrixGroup
end

function symplectic_doubling(A::MatrixElem)
    Ati = transpose(inv(A))
    return block_diagonal_matrix([A,Ati])
end

function SymplecticDoubling(G::MatrixGroup)
    GD = matrix_group([symplectic_doubling(matrix(g)) for g in gens(G)])
    return SymplecticDoubling(G, GD)
end

function SymplecticDoubling(W::Union{Gapjm.PRG, Gapjm.FiniteCoxeterGroup})
    return SymplecticDoubling(matrix_group(W))
end
