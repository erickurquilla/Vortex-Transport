#include <vector>
#include <iostream>

#include "Lagrangebasis.H"
#include "Profiling.H"

PROFILE_DECLARE("lagrange_basis_reference_space");
PROFILE_DECLARE("lagrange_basis_gradient_reference_space");
PROFILE_DECLARE("reference_to_physical_space");

// compute the lagrange basis function in reference space ( xi , eta )
std::vector<double> lagrange_basis_reference_space(const int& p, const std::vector<double>& coords_ref_spa){

    PROFILE_SCOPE("lagrange_basis_reference_space");

    // coordinates in reference space ( xi , eta )
    double xi, eta;
    xi = coords_ref_spa[0];
    eta = coords_ref_spa[1];
    
    // phi( xhi , eta ) for each interior node in the element 
    std::vector<double> phi( ( p + 1 ) *( p + 2 ) / 2 );

    switch(p) {
        case 0:
            phi[0] = 1.0;
            break;
        case 1:
            phi[0] = 1.0 - xi - eta;
            phi[1] = xi;
            phi[2] = eta;
            break;
        case 2:
            phi[0] = 1.0 - 3.0*xi - 3.0*eta + 2.0*xi*xi + 4.0*xi*eta + 2.0*eta*eta;
            phi[1] = 4.0*xi - 4.0*xi*xi - 4.0*xi*eta;
            phi[2] = -xi + 2.0*xi*xi;
            phi[3] = 4.0*eta - 4.0*xi*eta - 4.0*eta*eta;
            phi[4] = 4.0*xi*eta;
            phi[5] = -eta + 2.0*eta*eta;
            break;
        case 3:
            phi[0] = 1.0 - 11.0/2.0*xi - 11.0/2.0*eta + 9.0*xi*xi + 18.0*xi*eta + 9.0*eta*eta - 9.0/2.0*xi*xi*xi - 27.0/2.0*xi*xi*eta - 27.0/2.0*xi*eta*eta - 9.0/2.0*eta*eta*eta;
            phi[1] = 9.0*xi - 45.0/2.0*xi*xi - 45.0/2.0*xi*eta + 27.0/2.0*xi*xi*xi + 27.0*xi*xi*eta + 27.0/2.0*xi*eta*eta;
            phi[2] = -9.0/2.0*xi + 18.0*xi*xi + 9.0/2.0*xi*eta - 27.0/2.0*xi*xi*xi -27.0/2.0*xi*xi*eta;
            phi[3] = xi - 9.0/2.0*xi*xi + 9.0/2.0*xi*xi*xi;
            phi[4] = 9.0*eta - 45.0/2.0*xi*eta - 45.0/2.0*eta*eta + 27.0/2.0*xi*xi*eta + 27.0*xi*eta*eta + 27.0/2.0*eta*eta*eta;
            phi[5] = 27.0*xi*eta - 27.0*xi*xi*eta - 27.0*xi*eta*eta;
            phi[6] = -9.0/2.0*xi*eta + 27.0/2.0*xi*xi*eta;
            phi[7] = -9.0/2.0*eta + 9.0/2.0*xi*eta + 18.0*eta*eta - 27.0/2.0*xi*eta*eta - 27.0/2.0*eta*eta*eta;
            phi[8] = -9.0/2.0*xi*eta + 27.0/2.0*xi*eta*eta;
            phi[9] = eta - 9.0/2.0*eta*eta + 9.0/2.0*eta*eta*eta;
            break;
        default:
            printf("ERROR: Unsupported p order for lagrange_basis_reference_space");
            exit(EXIT_FAILURE);
            break;
    }

    // return the value of all lagrange basis function at position ( xi , eta ) in reference space
    // first index represent the node number, runs between 0 and ( p + 1 ) *( p + 2 ) / 2
    // the value contained in the lagran bases function evalueted in point ( xi , eta ) in reference space
    return phi;
}

// compute the gradiente of the lagrange basis function in reference space ( xi , eta )
std::vector<std::vector<double>> lagrange_basis_gradient_reference_space(const int& p, const std::vector<double>& coords_ref_spa){

    PROFILE_SCOPE("lagrange_basis_gradient_reference_space");

    // coordinates in reference space ( xi , eta )
    double xi, eta;
    xi = coords_ref_spa[0];
    eta = coords_ref_spa[1];

    // gradient of phi in position ( xhi , eta ) for each interior node in the element 
    // first index represent the node number, runs between 0 and ( p + 1 ) *( p + 2 ) / 2
    // second index represent physical space. 0: x and 1: y
    std::vector<std::vector<double>> gradient_phi( ( p + 1 ) *( p + 2 ) / 2 , std::vector<double>(2) );

    switch(p) {
        case 0:
            // x component
            gradient_phi[0][0] = 0.0;
            // y component
            gradient_phi[0][1] = 0.0;
            break;
        case 1:
            // x component
            gradient_phi[0][0] = -1.0;
            gradient_phi[1][0] = 1.0;
            gradient_phi[2][0] = 0.0;
            // y component
            gradient_phi[0][1] = -1.0;
            gradient_phi[1][1] = 0.0;
            gradient_phi[2][1] = 1.0;
            break;
        case 2:
            // x component
            gradient_phi[0][0] = -3.0 + 4.0*xi + 4.0*eta;
            gradient_phi[1][0] = 4.0 - 8.0*xi - 4.0*eta;
            gradient_phi[2][0] = -1.0 + 4.0*xi;
            gradient_phi[3][0] = -4.0 * eta;
            gradient_phi[4][0] = 4.0 * eta;
            gradient_phi[5][0] = 0.0;
            // y component
            gradient_phi[0][1] = -3.0 + 4.0*xi + 4.0*eta;
            gradient_phi[1][1] = -4.0*xi;
            gradient_phi[2][1] = 0.0;
            gradient_phi[3][1] = 4.0 - 4.0*xi - 8.0*eta;
            gradient_phi[4][1] = 4.0*xi;
            gradient_phi[5][1] = -1.0 + 4.0*eta;
            break;
        case 3:
            // x component    
            gradient_phi[0][0] = -11.0/2.0 + 18.0*xi + 18.0*eta - 27.0/2.0*xi*xi - 27.0*xi*eta - 27.0/2.0*eta*eta;
            gradient_phi[1][0] = 9.0 - 45.0*xi - 45.0/2.0*eta + 81.0/2.0*xi*xi + 54.0*xi*eta + 27.0/2.0*eta*eta;
            gradient_phi[2][0] = -9.0/2.0 + 36.0*xi + 9.0/2.0*eta - 81.0/2.0*xi*xi - 27.0*xi*eta;
            gradient_phi[3][0] = 1.0 - 9.0*xi + 27.0/2.0*xi*xi;
            gradient_phi[4][0] = -45.0/2.0*eta + 27.0*xi*eta + 27.0*eta*eta;
            gradient_phi[5][0] = 27.0*eta - 54.0*xi*eta - 27.0*eta*eta;
            gradient_phi[6][0] = -9.0/2.0*eta + 27.0*xi*eta;
            gradient_phi[7][0] = 9.0/2.0*eta - 27.0/2.0*eta*eta;
            gradient_phi[8][0] = -9.0/2.0*eta + 27.0/2.0*eta*eta;
            gradient_phi[9][0] = 0.0;
            // y component
            gradient_phi[0][1] = -11.0/2.0 + 18.0*xi + 18.0*eta - 27.0/2.0*xi*xi - 27.0*xi*eta - 27.0/2.0*eta*eta;
            gradient_phi[1][1] = -45.0/2.0*xi + 27.0*xi*xi + 27.0*xi*eta;
            gradient_phi[2][1] = 9.0/2.0*xi - 27.0/2.0*xi*xi;
            gradient_phi[3][1] = 0.0;
            gradient_phi[4][1] = 9.0 - 45.0/2.0*xi - 45.0*eta + 27.0/2.0*xi*xi + 54.0*xi*eta + 81.0/2.0*eta*eta;
            gradient_phi[5][1] = 27.0*xi - 27.0*xi*xi - 54.0*xi*eta;
            gradient_phi[6][1] = -9.0/2.0*xi + 27.0/2.0*xi*xi;
            gradient_phi[7][1] = -9.0/2.0 + 9.0/2.0*xi + 36.0*eta - 27.0*xi*eta - 81.0/2.0*eta*eta;
            gradient_phi[8][1] = -9.0/2.0*xi + 27.0*xi*eta;
            gradient_phi[9][1] = 1.0 - 9.0*eta + 27.0/2.0*eta*eta;
            break;
        default:
            printf("ERROR: Unsupported p order for lagrange_basis_gradient_reference_space");
            exit(EXIT_FAILURE);
            break;
    }

    // gradient of phi in position ( xhi , eta ) for each interior node in the element 
    // first index represent the node number, runs between 0 and ( p + 1 ) *( p + 2 ) / 2
    // second index represent physical space. 0: x and 1: y
    return gradient_phi;
}

// compute coordinates in physical space given ( xi , eta ) in reference space
std::vector<double> reference_to_physical_space(const std::vector<double>& coord_ref_spa, const std::vector<std::vector<double>>& vertex_phys_spa){

    PROFILE_SCOPE("reference_to_physical_space");

    // phi( xhi , eta )
    std::vector<double> phi = lagrange_basis_reference_space(1, coord_ref_spa);

    // interpolation from reference to physical space
    std::vector<double> r(2);

    r[0] = vertex_phys_spa[0][0] * phi[0] + vertex_phys_spa[1][0] * phi[1] + vertex_phys_spa[2][0] * phi[2]; // x coordinate
    r[1] = vertex_phys_spa[0][1] * phi[0] + vertex_phys_spa[1][1] * phi[1] + vertex_phys_spa[2][1] * phi[2]; // y coordinate

    return r;

}
