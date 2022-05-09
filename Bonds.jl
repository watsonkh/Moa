using LinearAlgebra
mutable struct Bond
    from::PD.Node
    to::PD.Node
    isBroken::Bool
end
Bond(from::PD.Node, to::PD.Node) = Bond(from, to, false)

function get_strain(bond::Bond)
    initial_bond_length = norm(bond.to.position - bond.from.position)
    deformed_bond_vector = (bond.to.position + bond.to.displacement) - (bond.from.position + bond.from.displacement)
    deformed_bond_length = norm(deformed_bond_vector)
    return (deformed_bond_length - initial_bond_length) / initial_bond_length
end

"Returns the force of the bond with the minimum material properties"
function get_force(bond::Bond)
    if bond.isBroken
        return zeros(3,)
    end
    # Duplicated code here from get_strain for efficiency
    initial_bond_length = norm(bond.to.position - bond.from.position)
    deformed_bond_vector = (bond.to.position + bond.to.displacement) - (bond.from.position + bond.from.displacement)
    deformed_bond_length = norm(deformed_bond_vector)
    strain = (deformed_bond_length - initial_bond_length) / initial_bond_length
    direction = deformed_bond_vector / deformed_bond_length

    return minimum([bond.from.material.bond_constant, bond.to.material.bond_constant]) * strain * direction * bond.to.volume * bond.from.volume
end

"Applies the bond's force to its from node (NOT THREAD SAFE)"
function apply_force(bond::Bond)
    @atomic bond.from.force += get_force(bond)
end

function should_break(bond::Bond)
    return get_strain(bond) > min(bond.from.material.critical_strain, bond.to.material.critical_strain)
end

function break!(bond::Bond)
    # println("SNAPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP")
    bond.isBroken = true
end

