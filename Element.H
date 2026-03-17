#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include <string>

#include "Meshgeneration.H"
#include "Parameters.H"

class Element {
public:

    // Default constructor
    Element();
 
    // Element default constructor with arguments
    Element(const int& ele_num, const mesh& mesh_info, const std::vector<std::vector<double>>& nods_ref_spa, const int& p_lagrange);

    // Element properties
    double time; // storage time in the dimulation for each element
    int p; // lagrange polinomial order
    int number; // element number
    int type; // element type 0 of squere angle is down and 1 if up
    int right_element; // elemnt to the right
    int left_element; // elemnt to the left
    int vertical_element; // elemnt in the vertical direction
    std::vector<std::vector<double>> vertices_coords_phys_space; // coordinates of the verctices of the element in physical space
    std::vector<std::vector<double>> nods_coords_refe_space; // coordinates of the interior nodes of the element in reference space
    std::vector<std::vector<double>> nods_coords_phys_space; // coordinates of the interior nodes of the element in physical space

    std::vector<std::vector<double>> hidrodynamics_vector_U; // containg the hydrodinamics vecto U for each interior node
    std::vector<std::vector<std::vector<double>>> hidrodynamics_vector_F; // containg the hydrodinamics vecto F for each interior node

    std::vector<std::vector<double>> jacobian; // jacobian between transformation from reference space to physical space d vec{x} / d vec{xi} = [ [ x2 - x1 , x3 - x1 ] , [ y2 - y1 , y3 - y1 ] ]
    double determinant_jacobian; // determinant of jacobian
    std::vector<std::vector<double>> inverse_jacobian; // jacobian between transformation from physical space to reference space d vec{xi} / d vec{x} = ( 1 / det( J ) ) * [ [ y3 - y1 , x1 - x3 ] , [ y1 - y2 , x2 - x1 ] ]
    double determinant_inverse_jacobian; // determinant of the inverse of the jacobian

    std::vector<std::vector<double>> inverse_mass_matrix_physical_space; // mass_ij = int in Omega phi_i phi_j dOmega . Size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
    std::vector<std::vector<std::vector<double>>> stiffness_matrix_physical_space; // S_ij = integral in Omega of ( Nabla phi_i ) phi_j dOmega . form :  hat{e}_x * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] + hat{e}_y * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ]. first index run between spacial components in physics space. 0: x and 1 y. second and third index run over matrix inidices of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] 
    
    std::vector<std::vector<double>> units_vectors_perpendicular_to_element_boundary; // contains the units vectors perperdicular to the elements boundary. the first index runs from 0 two 2 and represent the side number. side 1 is the one found going counterclockwise from the initial vextex (square angle vertex), in continuation side 2 and 3 going conter clockwise. the second index runs from 0 to 1, 0 is the x component of the unit vector 1 is the y components.

    std::vector<double> sides_lenght; // contains the element side lenghts, this array contains just three values. 0: side 1, 1: side 2, 2: side 3

    // Element functions
    void initialize_hydrodinamics(const int& ini_type, const std::vector<std::vector<double>>& gau_area_int); // initialize the hydronimics quantities U and F
    void write_data(const int& step_num); // write element data in output directory
    void build_mass_matrix_inverse(const std::vector<std::vector<double>>& inv_mass_matrix); // builds mass matrix inverse from referece space to physical space for each element
    void build_stiffness_matrix(const std::vector<std::vector<std::vector<double>>>& stiff_matrix); // builds stiffness matrix from referece space to physical space for each element
    void build_jacobians(); // compute jacobians to connetc referece space to physical space and viceversa for each element

};

#endif // ELEMENT_H