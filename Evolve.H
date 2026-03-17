#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

#include "Element.H"

class Evolve_element {
public:

    // Default constructor
    Evolve_element();
 
    // Element default constructor with arguments
    Evolve_element(Element* this_elem, Element* right_elem, Element* left_elem, Element* vertical_elem, const std::vector<std::vector<double>>& gau_int_l,  const int& p_parms);

    // Evolve_element properties
    int p; // lagrange polinomial order
    Element* element_this;     // this is the element to be evolved in time
    Element* element_right;    // this the element to the right to the element to be evolved in time
    Element* element_left;     // this the element to the left to the element to be evolved in time
    Element* element_vertical; // this the element in the vertical direction to the element to be evolved in time
    std::vector<std::vector<double>> gau_integ_line; // contains gauss quadrature information for line integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 1. 0 is xi coordinate. 1 is weight.

    // save the lagrange polinomial function evaluated in the quadrature point for the line integral for side 1, 2 and 3
    // this vector store the values of the lagrange polinomial in evaluated in the quadrature points
    // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
    // second item runs over the evaluation of the lagrange poliniam in the quadrature points
    // plus represent this element
    std::vector<std::vector<double>> plus_phi_in_quadrature_points_side_1; 
    std::vector<std::vector<double>> plus_phi_in_quadrature_points_side_2;
    std::vector<std::vector<double>> plus_phi_in_quadrature_points_side_3;
    // minus represent boundary element
    std::vector<std::vector<double>> minus_phi_in_quadrature_points_side_1; 
    std::vector<std::vector<double>> minus_phi_in_quadrature_points_side_2;
    std::vector<std::vector<double>> minus_phi_in_quadrature_points_side_3;

    // the following vectors store the interpolation of the hidrodynamic vector u from the interior nodes of the elements to the quadrature points for the line integrals
    // first index runs over the quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
    // plus for current element
    std::vector<std::vector<double>> U_plus_side_1; 
    std::vector<std::vector<double>> U_plus_side_2;
    std::vector<std::vector<double>> U_plus_side_3;
    // minus for elements in the boundaries
    std::vector<std::vector<double>> U_minus_side_1;
    std::vector<std::vector<double>> U_minus_side_2;
    std::vector<std::vector<double>> U_minus_side_3;

    // contains the numerical flux (roe flux) at the quadrature points for sides 1, 2 and 3
    // first index runs over the gauss quadrature points number, second index runs over 0 and 3 for the hidrodynamics quantities
    std::vector<std::vector<double>> numerical_flux_side_1;
    std::vector<std::vector<double>> numerical_flux_side_2;
    std::vector<std::vector<double>> numerical_flux_side_3;

    // Discontinuos Galerkin method variables

    // DG vector that results from the integration of the numerical flux ( integral phi_i hat_{F} dl ). 
    // First index runs over interior nodes. 
    // Second index runs between 0 and 3 and represend hidrodynamics variables.
    std::vector<std::vector<double>> DG_numerical_flux_integration; 
    // Stiffness vector: DG vector that results from area integral over element of nabla phi_i dot F dOmega.
    // First index runs over interior nodes. 
    // Second index runs between 0 and 3 and represend hidrodynamics variables.
    std::vector<std::vector<double>> DG_stiffness_vector; 
    // residual vector: DG vector that results from stiffness vector minus vector result of the numerical flux integration
    // First index runs over interior nodes. 
    // Second index runs between 0 and 3 and represend hidrodynamics variables.
    std::vector<std::vector<double>> DG_residual_vector; 
    // Time_derivative of state vector U(t): DG vector that results from dot{U}_{i} = M^{-1}_{ij} R_{j} ---> mass matrix inverse times the residual vector.
    // First index runs over interior nodes. 
    // Second index runs between 0 and 3 and represend hidrodynamics variables.
    std::vector<std::vector<double>> DG_time_derivative_U; 

    // Evolve_element functions
    void evaluate_basis_in_quadrature_poits(); // evaluate lagrande polinomials in the gauss quadrature at the element boundaries
    void compute_U_plus_minus(); // this function compute the hidrodynamic vector U on the element boundaries, side 1, 2 and 3. 
    void compute_numerical_flux(); // this function compute the numerical flux on the element boundaries, side 1, 2 and 3. 
    void integrate_numerical_flux(); // this function create the DG vector that results from the integration of the numerical flux ( integral phi_i hat_{F} dl ), this function fill the vector DG_numerical_flux_integration.
    void compute_stiffness_vector(); // compute stiffness vector (area integral over element of nabla phi_i dot F dOmega)
    void compute_residual_vector(); // compute residial vector stiffness vector minus vector result of the numerical flux integration
    void compute_time_derivative_U(); // compute DG_time_derivative_U that is time_derivative of state vector U(t): DG vector that results from dot{U}_{i} = M^{-1}_{ij} R_{j} ---> mass matrix inverse times the residual vector.

};

#endif // EVOLVE_H