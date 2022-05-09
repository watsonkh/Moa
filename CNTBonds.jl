using LinearAlgebra


import CSV
global strainForce = CSV.File("strain-force.csv",header=false)


mutable struct CNTBond
    from::PD.Node
    to::PD.Node
    isBroken::Bool
end

CNTBond(from::PD.Node, to::PD.Node) = CNTBond(from, to, false)

function get_strain(bond::CNTBond)
    initial_bond_length = norm(bond.to.position - bond.from.position)
    deformed_bond_vector = (bond.to.position + bond.to.displacement) - (bond.from.position + bond.from.displacement)
    deformed_bond_length = norm(deformed_bond_vector)
    return (deformed_bond_length - initial_bond_length) / initial_bond_length
end

"Returns the force of the bond with the minimum material properties"
function get_force(bond::CNTBond)
    if bond.isBroken
        return zeros(3,)
    end
    # Duplicated code here from get_strain for efficiency
    initial_bond_length = norm(bond.to.position - bond.from.position)
    deformed_bond_vector = (bond.to.position + bond.to.displacement) - (bond.from.position + bond.from.displacement)
    deformed_bond_length = norm(deformed_bond_vector)
    strain = (deformed_bond_length - initial_bond_length) / initial_bond_length
    direction = deformed_bond_vector / deformed_bond_length
    
    # Linearly interpolate and extrapolate if out of bounds
    custom_force_coeff = 1
    if strain < strainForce[1][1]
        custom_force_coeff = (strainForce[2][2] - strainForce[1][2]) * (strain - strainForce[1][1]) / (strainForce[2][1] - strainForce[1][1]) + strainForce[1][2]
    elseif strain > strainForce[length(strainForce)][1]
        custom_force = (strainForce[length(strainForce)][2] - strainForce[length(strainForce)-1][2]) * (strain - strainForce[length(strainForce)-1][1]) / (strainForce[length(strainForce)][1] - strainForce[length(strainForce)-1][1]) + strainForce[length(strainForce)-1][2]
    else
        for i in 2:length(strainForce)
            if strain < strainForce[i][1]
                custom_force = (strainForce[i][2] - strainForce[i-1][2]) * (strain - strainForce[i-1][1]) / (strainForce[i][1] - strainForce[i-1][1]) + strainForce[i-1][2]
                break
            end
        end
    end
    # println(strain, ", ", custom_force)
    return custom_force / 6. * direction
end

"Applies the bond's force to its from node (NOT THREAD SAFE)"
function apply_force(bond::CNTBond)
    @atomic bond.from.force += get_force(bond)
end

function should_break(bond::CNTBond)
    return get_strain(bond) > min(bond.from.material.critical_strain, bond.to.material.critical_strain)
end

function break!(bond::CNTBond)
    # println("SNAPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP")
    bond.isBroken = true
end
