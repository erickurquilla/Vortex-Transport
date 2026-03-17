#include <iostream>
#include <vector>

#include "Parameters.H"
#include "Utilities.H"
#include "Meshgeneration.H"
#include "Element.H"
#include "Quadraturerule.H"
#include "Preevolve.H"
#include "Evolve.H"
#include "Timestepping.H"

int main(int argc, char* argv[]) {
    
    // read simulation paramaters
    parameters parms = read_input_files(argc, argv);

    // generate mesh
    mesh simulation_mesh = generate_mesh(parms);
    
    // generate gauss quadrature for line and area integral
    std::vector<std::vector<double>> gauss_integral_line = gauss_line_integral(parms.integration_order);
    std::vector<std::vector<double>> gauss_integral_area = gauss_area_integral(parms.integration_order);

    // generate elements interior nodes in reference space
    std::vector<std::vector<double>> nodes_reference_space = generate_nodes_reference_space(parms);

    // compute the inverse of the mass matrix in reference space
    std::vector<std::vector<double>> inverse_mass_matrix = inverse_mass_matrix_reference_space(parms.p, gauss_integral_area);

    // compute stiffness matrix in reference space
    std::vector<std::vector<std::vector<double>>> stiffness_matrix = sitffness_matrix_reference_space(parms.p, gauss_integral_area);

    // Define an array of Elements members of the class Elements
    Element* elements = new Element[  2 * parms.num_element_in_x * parms.num_element_in_y  ];

    // Create the output directory
    clean_create_directory("output");
    clean_create_directory("output/step_0");

    // Initialize elements of the array
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y ; ++i) {
        elements[i] = Element(i, simulation_mesh, nodes_reference_space, parms.p); // Initilize element objects and compute interior nodes coordinate in physical space
        elements[i].build_jacobians(); // compute jacobians to connetc referece space to physical space and viceversa for each element
        elements[i].build_mass_matrix_inverse(inverse_mass_matrix); // Build mass matrix
        elements[i].build_stiffness_matrix(stiffness_matrix); // builds stiffness matrix from referece space to physical space for each element
        elements[i].initialize_hydrodinamics(parms.U_initialization_type, gauss_integral_area); //Initialize hidrodynamic quanities
        elements[i].write_data(0); //write data
    }

    // Define an array of Evolve_elements members of the class Evolve_element
    Evolve_element* evolve_elements = new Evolve_element[  2 * parms.num_element_in_x * parms.num_element_in_y  ];

    // Initialize Evolve_element objects in the array evolve_elements
    for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y ; ++i) {
        // Build evolve_elements
        evolve_elements[i] = Evolve_element(&elements[i],&elements[elements[i].right_element],&elements[elements[i].left_element],&elements[elements[i].vertical_element], gauss_integral_line, parms.p);
        evolve_elements[i].evaluate_basis_in_quadrature_poits();
    }

    // time stepping loop
    for (int a = 1; a < parms.number_time_steps + 1; ++a) {
        // 0 : forward euler
        if ( parms.stepping_method == 0 ){
            // called forward_euler function in Timestepping.cpp script
            forward_euler(elements, evolve_elements, parms.time_step, 2 * parms.num_element_in_x * parms.num_element_in_y);
        // 1 : rk4
        }else if ( parms.stepping_method == 1 ){
            // called rk4 function in Timestepping.cpp script
            rk4(elements, evolve_elements, parms.time_step, 2 * parms.num_element_in_x * parms.num_element_in_y);
        }else{
            printf("ERROR: Unsupported time stepping method\n0 : forward euler\n1 : rk4\n");
            exit(EXIT_FAILURE);
        }

        // write data
        if (a % parms.write_every_steps == 0) {
            // Create the output/step_a directory
            clean_create_directory("output/step_" + std::to_string(a));        
            std::cout << "Step : " << a << std::endl;
            // loop over elements
            for (int i = 0; i < 2 * parms.num_element_in_x * parms.num_element_in_y ; ++i) {
                    // write data
                    elements[i].write_data(a); 
            }
        }
    }

    return 0;
}