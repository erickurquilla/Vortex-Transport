#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Lagrangebasis.H"
#include "Utilities.H"

Element::Element() {}

Element::Element(const int& ele_num, const mesh& mesh_info, const std::vector<std::vector<double>>& nods_ref_spa, const int& p_lagrange):
    
    // Initialize element properties    
    time(0.0), // Simulation initial time
    p(p_lagrange), // lagrange polinomial order
    number(ele_num), // Element number
    type(mesh_info.element_type[ele_num]), // Element type 0 of squere angle is down and 1 if up
    right_element(mesh_info.elements_at_boundary[ele_num][0]), // Elemnt to the right
    left_element(mesh_info.elements_at_boundary[ele_num][1]), // Elemnt to the left
    vertical_element(mesh_info.elements_at_boundary[ele_num][2]), // Elemnt in the vertical direction
    vertices_coords_phys_space(mesh_info.element_coordinates[ele_num]), // Coordinates of the verctices of the element in physical space. First index represent vertices of the element (run between 0 and 2): 0 is the vertice at the square angle, 1 and 2 are the other vertices going counter clock wise. Second index represent x:0 and y:1 position
    nods_coords_refe_space(nods_ref_spa), // Coordinates of the interior nodes of the element in reference space. The first index in the node number, the second inxed runs between 0 and 1. 0: xi and 1: eta.
    nods_coords_phys_space((this->p + 1) * (this->p + 2) / 2, std::vector<double>(2)), // Coordinates of the interior nodes of the element in physical space. The first index is the node number, the second inxed runs between 0 and 1. 0: x and 1: y.
    hidrodynamics_vector_U((this->p + 1) * (this->p + 2) / 2, std::vector<double>(4)), // Hidrodynamics: vector u for each interior node. The first index in the node number, the second inxed runs between 0 and 4 for hidrodynamics quantities. 0: rho, 1: rho u, 2: rho v, 3: rho E.
    hidrodynamics_vector_F((this->p + 1) * (this->p + 2) / 2, std::vector<std::vector<double>>(2, std::vector<double>(4))), // Hidrodynamics: vector f for each interior node. The first index in the node number. The second index runs between 0 and 1 and represente physical space 0:x and 1: y. The third inxed runs between 0 and 4 for hidrogynamics quantities.
    jacobian(2,std::vector<double>(2)), // jacobian between transformation from reference space to physical space d vec{x} / d vec{xi} = [ [ x2 - x1 , x3 - x1 ] , [ y2 - y1 , y3 - y1 ] ]
    inverse_jacobian(2,std::vector<double>(2)), // jacobian between transformation from physical space to reference space d vec{xi} / d vec{x} = ( 1 / det( J ) ) * [ [ y3 - y1 , x1 - x3 ] , [ y1 - y2 , x2 - x1 ] ]
    inverse_mass_matrix_physical_space((this->p + 1) * (this->p + 2) / 2 , std::vector<double>((this->p + 1) * (this->p + 2) / 2)), // mass_ij = int in Omega phi_i phi_j dOmega . Size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
    stiffness_matrix_physical_space(2, std::vector<std::vector<double>>( (this->p + 1) * (this->p + 2) / 2 , std::vector<double>( (this->p + 1) * (this->p + 2) / 2 ))), // S_ij = integral in Omega of ( Nabla phi_i ) phi_j dOmega . form :  hat{e}_x * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] + hat{e}_y * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ]. first index run between spacial components in physics space. 0: x and 1 y. second and third index run over matrix inidices of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] 
    units_vectors_perpendicular_to_element_boundary(3,std::vector<double>(2)), // contains the units vectors perperdicular to the elements boundary. the first index runs from 0 two 2 and represent the side number. side 1 is the one found going counterclockwise from the initial vextex (square angle vertex), in continuation side 2 and 3 going conter clockwise. the second index runs from 0 to 1, 0 is the x component of the unit vector 1 is the y components.
    sides_lenght(3) // contains the element side lenghts, this array contains just three values. 0: side 1, 1: side 2, 2: side 3

    {
    
    // compute interior nodes coordinate in physical space
    for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
        this->nods_coords_phys_space[i] = reference_to_physical_space(this->nods_coords_refe_space[i], this->vertices_coords_phys_space);
    }
}

// compute jacobians to connetc referece space to physical space and viceversa for each element
void Element::build_jacobians(){

    // compute jacobian between transformation from reference space to physical space d vec{x} / d vec{xi} = [ [ x2 - x1 , x3 - x1 ] , [ y2 - y1 , y3 - y1 ] ] 
    this->jacobian[0][0] = this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0];
    this->jacobian[0][1] = this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[0][0];
    this->jacobian[1][0] = this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1];
    this->jacobian[1][1] = this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[0][1];

    // compute determinant of jacobian
    this->determinant_jacobian = this->jacobian[0][0] * this->jacobian[1][1] - this->jacobian[1][0] * this->jacobian[0][1];

    // compute jacobian between transformation from physical space to reference space d vec{xi} / d vec{x} = ( 1 / det( J ) ) * [ [ y3 - y1 , x1 - x3 ] , [ y1 - y2 , x2 - x1 ] ]
    this->inverse_jacobian[0][0] = ( 1 / this->determinant_jacobian ) * ( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[0][1] );
    this->inverse_jacobian[0][1] = ( 1 / this->determinant_jacobian ) * ( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] );
    this->inverse_jacobian[1][0] = ( 1 / this->determinant_jacobian ) * ( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[1][1] );
    this->inverse_jacobian[1][1] = ( 1 / this->determinant_jacobian ) * ( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] );

    // compute determinant of inverse jacobian
    this->determinant_inverse_jacobian = this->inverse_jacobian[0][0] * this->inverse_jacobian[1][1] - this->inverse_jacobian[1][0] * this->inverse_jacobian[0][1];

    // compute unit vectors perpendicular to the element edges
    //                                                      (            y2                     -              y1                 ) / (       (             x2                   -           x1                     )^2   +    (             y2                   -           y1                     )^2   )^1/2             
    this->units_vectors_perpendicular_to_element_boundary[0][0] = ( this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1] ) / pow( pow( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] , 2 ) + pow( this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1] , 2 ) , 0.5 );
    this->units_vectors_perpendicular_to_element_boundary[1][0] = ( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[1][1] ) / pow( pow( this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[1][0] , 2 ) + pow( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[1][1] , 2 ) , 0.5 );
    this->units_vectors_perpendicular_to_element_boundary[2][0] = ( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[2][1] ) / pow( pow( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] , 2 ) + pow( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[2][1] , 2 ) , 0.5 );

    //                                                      -1 * (            x2                     -              x1                 ) / (       (             x2                   -           x1                     )^2   +    (             y2                   -           y1                     )^2   )^1/2             
    this->units_vectors_perpendicular_to_element_boundary[0][1] = -1 * ( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] ) / pow( pow( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] , 2 ) + pow( this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1] , 2 ) , 0.5 );
    this->units_vectors_perpendicular_to_element_boundary[1][1] = -1 * ( this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[1][0] ) / pow( pow( this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[1][0] , 2 ) + pow( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[1][1] , 2 ) , 0.5 );
    this->units_vectors_perpendicular_to_element_boundary[2][1] = -1 * ( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] ) / pow( pow( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] , 2 ) + pow( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[2][1] , 2 ) , 0.5 );

    // compute the element side lenghts
    //                   (     (           x2                    -             x1                   )^2   +    (             y2                   -              y1                  )^2   )^0.5      
    this->sides_lenght[0] = pow( pow( this->vertices_coords_phys_space[1][0] - this->vertices_coords_phys_space[0][0] , 2 ) + pow( this->vertices_coords_phys_space[1][1] - this->vertices_coords_phys_space[0][1] , 2 ) , 0.5 ); // side 1
    this->sides_lenght[1] = pow( pow( this->vertices_coords_phys_space[2][0] - this->vertices_coords_phys_space[1][0] , 2 ) + pow( this->vertices_coords_phys_space[2][1] - this->vertices_coords_phys_space[1][1] , 2 ) , 0.5 ); // side 2
    this->sides_lenght[2] = pow( pow( this->vertices_coords_phys_space[0][0] - this->vertices_coords_phys_space[2][0] , 2 ) + pow( this->vertices_coords_phys_space[0][1] - this->vertices_coords_phys_space[2][1] , 2 ) , 0.5 ); // side 3

}

// builds mass matrix inverse from referece space to physical space for each element
void Element::build_mass_matrix_inverse(const std::vector<std::vector<double>>& inv_mass_matrix){
    // compute inverse mass matrix in physical space
    for (int i = 0; i < (this->p + 1) * (this->p + 2) / 2; ++i) {
        for (int j = 0; j < (this->p + 1) * (this->p + 2) / 2; ++j) {
            this->inverse_mass_matrix_physical_space[i][j] = inv_mass_matrix[i][j] / this->determinant_jacobian;     
        }
    }
}

// builds stiffness matrix from referece space to physical space for each element
void Element::build_stiffness_matrix(const std::vector<std::vector<std::vector<double>>>& stiff_matrix){
    // compute stiffness matrix in physical space 
    for (int i = 0; i < (this->p + 1) * (this->p + 2) / 2; ++i) {
        for (int j = 0; j < (this->p + 1) * (this->p + 2) / 2; ++j) {
            this->stiffness_matrix_physical_space[0][i][j] = this->determinant_jacobian * ( stiff_matrix[0][i][j] * this->inverse_jacobian[0][0] + stiff_matrix[1][i][j] * this->inverse_jacobian[1][0] );     
            this->stiffness_matrix_physical_space[1][i][j] = this->determinant_jacobian * ( stiff_matrix[0][i][j] * this->inverse_jacobian[0][1] + stiff_matrix[1][i][j] * this->inverse_jacobian[1][1] );
        }
    }
}

// initialize the hydronimics quantities U and F
void Element::initialize_hydrodinamics(const int& ini_type, const std::vector<std::vector<double>>& gau_area_int){
    // initialization type of hidrodynamics state U
    // 0 : direct interpolation
    if ( ini_type == 0 ){
        // hidrodynamic quantities 
        double rho, u, v, p, E, H;
        // position (x,y) and time
        double x, y, t;
        // hidrodynamic constant
        double rho_infty, rc, epsilon, gamma, M_infty, p_infty, U_infty, V_infty, x0, y0;
        // initial conditions functions
        double f0, f1, f2;

        rho_infty = 1.0;
        rc = 1.0;
        epsilon = 0.3;
        gamma = 1.4;
        M_infty = 0.5;
        p_infty = 20.0 / 7.0;
        U_infty = 1.0 / pow( 2.0 , 0.5);
        V_infty = 1.0 / pow( 2.0 , 0.5);
        x0 = 0.0;
        y0 = 0.0;

        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {

            x = this->nods_coords_phys_space[i][0]; // x position of node i
            y = this->nods_coords_phys_space[i][1]; // y position of node i
            t = this->time; // initial time

            f0 = 1.0 - ( pow ( x - x0 - U_infty * t , 2.0 ) + pow ( y - y0 - V_infty * t , 2.0 ) ) / pow ( rc , 2.0 );
            f1 = 1.0 - pow ( epsilon , 2.0 ) * ( gamma -1 ) * pow ( M_infty , 2.0) * exp( f0 ) / ( 8.0 * pow ( M_PI , 2.0 ) );
            f2 = epsilon * ( pow( U_infty , 2.0) + pow( V_infty , 2.0) ) * exp( f0 / 2.0 ) / ( 2.0 * M_PI * rc );

            rho = rho_infty * pow( f1 , 1.0 / ( gamma - 1.0 ) ); // density
            u   = U_infty - f2 * ( y - y0 - V_infty * t ); // horizontal velocity
            v   = V_infty + f2 * ( x - x0 - U_infty * t ); // vertical velocity
            p   = p_infty * pow( f1 , gamma / ( gamma - 1.0 ) ); // pressure
            E   = p / ( rho * ( gamma - 1.0 ) ) + ( pow( u , 2.0) + pow( v , 2.0) ) / 2.0; // Energy
            H   = E + p / rho; // Entalpy

            // initialize the hidrodinamic vector u  
            this->hidrodynamics_vector_U[i][0] = rho;
            this->hidrodynamics_vector_U[i][1] = rho * u;
            this->hidrodynamics_vector_U[i][2] = rho * v;
            this->hidrodynamics_vector_U[i][3] = rho * E;

            // initialize the hidrodinamic vector f, x component  
            this->hidrodynamics_vector_F[i][0][0] = rho * u;
            this->hidrodynamics_vector_F[i][0][1] = rho * pow( u , 2 ) + p;
            this->hidrodynamics_vector_F[i][0][2] = rho * u * v;
            this->hidrodynamics_vector_F[i][0][3] = rho * u * H;

            // initialize the hidrodinamic vector f, y component  
            this->hidrodynamics_vector_F[i][1][0] = rho * v;
            this->hidrodynamics_vector_F[i][1][1] = rho * u * v;
            this->hidrodynamics_vector_F[i][1][2] = rho * pow( v , 2 ) + p;
            this->hidrodynamics_vector_F[i][1][3] = rho * v * H;
        }
    }else if( ini_type == 1 ){ // 1 : least squere projection

        // number of quadrature points
        int number_quadrature_points = gau_area_int.size(); 

        // this vector store the values of the lagrange polinomial in evaluated in the quadrature points
        // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
        // second item runs over the evaluation of the lagrange poliniam in the quadrature points
        std::vector<std::vector<double>> phi_in_quad_points((this->p + 1) * (this->p + 2) / 2, std::vector<double>(number_quadrature_points));

        // this vector store the position of the quadrature points in physical space
        // first index runs over quadrature points
        // second item runs between 0 and 1. 0 is x position and 1 is the y position
        std::vector<std::vector<double>> r_phys_space_quad_points(number_quadrature_points, std::vector<double>(2));

        // position in reference space of quadrature point | temporal variable    
        std::vector<double> xi_eta_gauss = { 0.0 , 0.0 };
        // phi value in point xi_eta_gauss | temporal variable    
        std::vector<double> phi_in_xi_eta_gauss( (this->p + 1) * (this->p + 2) / 2 );
        

        // evaluate the lagrange polinomial in the quadrature points
        for (int i = 0; i < number_quadrature_points; ++i) {
            xi_eta_gauss[0] = gau_area_int[i][0]; // xi
            xi_eta_gauss[1] = gau_area_int[i][1]; // eta
            // evaluate phi value in point xi_eta_gauss
            phi_in_xi_eta_gauss = lagrange_basis_reference_space( this->p , xi_eta_gauss ); 
            // save the values of the lagrange polinomial in the quadrature points
            for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
                phi_in_quad_points[j][i] = phi_in_xi_eta_gauss[j];
            }
 
            // convert point xi_eta_gauss to physical space
            r_phys_space_quad_points[i] = reference_to_physical_space(xi_eta_gauss, this->vertices_coords_phys_space);

        }

        // DG indices vector u exact
        // first index runs over quadrature points
        // second item runs between 0 and 3 for hidrodynamic index
        std::vector<std::vector<double>> DG_u_exact(number_quadrature_points, std::vector<double>(4));

        // hidrodynamic quantities 
        double rho, u, v, p, E, H;
        // position (x,y) and time
        double x, y, t;
        // hidrodynamic constant
        double rho_infty, rc, epsilon, gamma, M_infty, p_infty, U_infty, V_infty, x0, y0;
        // initial conditions functions
        double f0, f1, f2;

        rho_infty = 1.0;
        rc = 1.0;
        epsilon = 0.3;
        gamma = 1.4;
        M_infty = 0.5;
        p_infty = 20.0 / 7.0;
        U_infty = 1.0 / pow( 2.0 , 0.5);
        V_infty = 1.0 / pow( 2.0 , 0.5);
        x0 = 0.0;
        y0 = 0.0;

        for (int i = 0; i < number_quadrature_points ; ++i) {

            x = r_phys_space_quad_points[i][0]; // x position of node i
            y = r_phys_space_quad_points[i][1]; // y position of node i
            t = time; // initial time

            f0 = 1.0 - ( pow ( x - x0 - U_infty * t , 2.0 ) + pow ( y - y0 - V_infty * t , 2.0 ) ) / pow ( rc , 2.0 );
            f1 = 1.0 - pow ( epsilon , 2.0 ) * ( gamma -1 ) * pow ( M_infty , 2.0) * exp( f0 ) / ( 8.0 * pow ( M_PI , 2.0 ) );
            f2 = epsilon * ( pow( U_infty , 2.0) + pow( V_infty , 2.0) ) * exp( f0 / 2.0 ) / ( 2.0 * M_PI * rc );

            rho = rho_infty * pow( f1 , 1.0 / ( gamma - 1.0 ) ); // density
            u   = U_infty - f2 * ( y - y0 - V_infty * t ); // horizontal velocity
            v   = V_infty + f2 * ( x - x0 - U_infty * t ); // vertical velocity
            p   = p_infty * pow( f1 , gamma / ( gamma - 1.0 ) ); // pressure
            E   = p / ( rho * ( gamma - 1.0 ) ) + ( pow( u , 2.0) + pow( v , 2.0) ) / 2.0; // Energy
            H   = E + p / rho; // Entalpy

            // initialize the hidrodinamic vector u  
            DG_u_exact[i][0] = rho;
            DG_u_exact[i][1] = rho * u;
            DG_u_exact[i][2] = rho * v;
            DG_u_exact[i][3] = rho * E;

        }

        // DG indices vector b
        // First index runs over interior nodes. 
        // Second index runs between 0 and 3 and represend hidrodynamics variables.
        std::vector<std::vector<double>> DG_b((this->p + 1) * (this->p + 2) / 2, std::vector<double>(4));

        // loop over all interior nodes
        for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {

            // initialize vector b at zero
            DG_b[i] = { 0.0 , 0.0, 0.0 , 0.0 };
            
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                // loop over quadrature points
                for (int j = 0; j < number_quadrature_points; ++j) {
                    DG_b[i][k] += this->determinant_jacobian * phi_in_quad_points[i][j] * DG_u_exact[j][k] * gau_area_int[j][2];
                    // DG_b is the DG vector that results from the integration of the u exect ( integral phi_i u_exact dOmega ). 
                    // First index runs over interior nodes. 
                    // Second index runs between 0 and 3 and represend hidrodynamics variables.
                }
            }
        }

        // compute hidrodynamic vector U
        // loop over all the interior nodes of this element
        for (int i = 0; i < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++i) {

            // initialize hidrodynamics_vector_U[i][j] values to zero
            this->hidrodynamics_vector_U[i] = { 0.0 , 0.0 , 0.0 , 0.0 }; 

            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                // loop over all the interior nodes of this element
                for (int j = 0; j < ( this->p + 1 ) * ( this->p + 2 ) / 2; ++j) {
                    //             U_i                  = sum_j              M^{-1}_{ij}                  *  b_{j}  
                    this->hidrodynamics_vector_U[i][k] += this->inverse_mass_matrix_physical_space[i][j] * DG_b[j][k]; 
                }
            }      

            rho = this->hidrodynamics_vector_U[i][0]; // density
            u   = this->hidrodynamics_vector_U[i][1] / rho; // horizontal velocity
            v   = this->hidrodynamics_vector_U[i][2] / rho; // vertical velocity
            E   = this->hidrodynamics_vector_U[i][3] / rho; // energy
            p   = rho * ( gamma - 1 ) * ( E - ( pow( u , 2) + pow( v , 2) ) / 2 ); // pressure
            H   = E + p / rho; // Entalpy

            // initialize the hidrodinamic vector f, x component  
            this->hidrodynamics_vector_F[i][0][0] = rho * u;
            this->hidrodynamics_vector_F[i][0][1] = rho * pow( u , 2 ) + p;
            this->hidrodynamics_vector_F[i][0][2] = rho * u * v;
            this->hidrodynamics_vector_F[i][0][3] = rho * u * H;

            // initialize the hidrodinamic vector f, y component  
            this->hidrodynamics_vector_F[i][1][0] = rho * v;
            this->hidrodynamics_vector_F[i][1][1] = rho * u * v;
            this->hidrodynamics_vector_F[i][1][2] = rho * pow( v , 2 ) + p;
            this->hidrodynamics_vector_F[i][1][3] = rho * v * H;

        }

    }else {
        printf("ERROR: Unsupported initialization type\n0 : direct interpolation\n1 : least squere projection\n");
        exit(EXIT_FAILURE);
    }  


}

// write element data in output directory
void Element::write_data(const int& step_num){

    // prepare data to be saved
    std::vector<std::string> lines( 1 + ( this->p + 1 ) *( this->p + 2 ) / 2 );
    lines[0]="node_number time x y u0 u1 u2 u3 fx0 fx1 fx2 fx3 fy0 fy1 fy2 fy3";

    for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {

        lines[i+1]=std::to_string(i)+" "; // node number
        lines[i+1]+=std::to_string(this->time)+" "; // time
        lines[i+1]+=std::to_string(this->nods_coords_phys_space[i][0])+" "; // x
        lines[i+1]+=std::to_string(this->nods_coords_phys_space[i][1])+" "; // y 
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_U[i][0])+" "; // u0
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_U[i][1])+" "; // u1
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_U[i][2])+" "; // u2
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_U[i][3])+" "; // u3
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_F[i][0][0])+" "; // fx0
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_F[i][0][1])+" "; // fx1
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_F[i][0][2])+" "; // fx2
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_F[i][0][3])+" "; // fx3
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_F[i][1][0])+" "; // fy0
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_F[i][1][1])+" "; // fy1
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_F[i][1][2])+" "; // fy2
        lines[i+1]+=std::to_string(this->hidrodynamics_vector_F[i][1][3])+" "; // fy3

    }

    writeToFile("output/step_" + std::to_string( step_num ) + "/element_" + std::to_string( this->number) + ".txt", lines);

    // if step is zero it writes the jacobians and their determinants, inverse mass and stiffness matrix
    if ( step_num == 0 ){

        // prepare jacobians and their determinants, inverse mass and stiffness matrix to be saved

        std::vector<std::string> lines_( 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 + 3 );
        
        lines_[0]="jacobian";
        lines_[1]="[ [ " + std::to_string( this->jacobian[0][0] ) + " , " + std::to_string( this->jacobian[0][1] ) + " ] , ";
        lines_[2]="  [ " + std::to_string( this->jacobian[1][0] ) + " , " + std::to_string( this->jacobian[1][1] ) + " ] ]";
        lines_[3]="determinant of jacobian";
        lines_[4]=std::to_string( this->determinant_jacobian );
        lines_[5]="inverse jacobian";
        lines_[6]="[ [ " + std::to_string( this->inverse_jacobian[0][0] ) + " , " + std::to_string( this->inverse_jacobian[0][1] ) + " ] , ";
        lines_[7]="  [ " + std::to_string( this->inverse_jacobian[1][0] ) + " , " + std::to_string( this->inverse_jacobian[1][1] ) + " ] ]";
        lines_[8]="determinant of inverse jacobian";
        lines_[9]=std::to_string( this->determinant_inverse_jacobian );
        lines_[10]="inverse mass matrix";
        lines_[11]="[ ";
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
            lines_[11+i]+="[ " + std::to_string(this->inverse_mass_matrix_physical_space[i][0]);
            for (int j = 1; j < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++j) {
                lines_[11+i]+=", " + std::to_string(this->inverse_mass_matrix_physical_space[i][j]);
            }
            lines_[11+i]+=" ] ,";
        }
        lines_[11 + ( this->p + 1 ) *( this->p + 2 ) / 2 -1 ].erase(lines_[11 + ( this->p + 1 ) *( this->p + 2 ) / 2 -1 ].size() - 1);
        lines_[11 + ( this->p + 1 ) *( this->p + 2 ) / 2 -1 ] += "]";
        lines_[11 + ( this->p + 1 ) *( this->p + 2 ) / 2]="stiffness matrix : x component";
        lines_[12 + ( this->p + 1 ) *( this->p + 2 ) / 2]="[ ";
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
            lines_[12 + ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+="[ " + std::to_string(this->stiffness_matrix_physical_space[0][i][0]);
            for (int j = 1; j < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++j) {
                lines_[12 + ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+=", " + std::to_string(this->stiffness_matrix_physical_space[0][i][j]);
            }
            lines_[12 + ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+=" ] ,";
        }
        lines_[12 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ].erase(lines_[ 12 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ].size() - 1);
        lines_[12 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] += "]";
        lines_[12 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2]="stiffness matrix : y component";
        lines_[13 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2]="[ ";
        for (int i = 0; i < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++i) {
            lines_[13 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+="[ " + std::to_string(this->stiffness_matrix_physical_space[1][i][0]);
            for (int j = 1; j < ( this->p + 1 ) *( this->p + 2 ) / 2 ; ++j) {
                lines_[13 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+=", " + std::to_string(this->stiffness_matrix_physical_space[1][i][j]);
            }
            lines_[13 + 2 * ( this->p + 1 ) *( this->p + 2 ) / 2 + i]+=" ] ,";
        }
        lines_[13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ].erase(lines_[ 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ].size() - 1);
        lines_[13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] += "]";

        lines_[1 + 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] = "Unit vector perpendicular to edge 1 : [ " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[0][0] ) + " , " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[0][1] ) + " ]";
        lines_[2 + 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] = "Unit vector perpendicular to edge 2 : [ " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[1][0] ) + " , " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[1][1] ) + " ]";
        lines_[3 + 13 + 3 * ( this->p + 1 ) *( this->p + 2 ) / 2 - 1 ] = "Unit vector perpendicular to edge 3 : [ " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[2][0] ) + " , " + std::to_string( this->units_vectors_perpendicular_to_element_boundary[2][1] ) + " ]";

        writeToFile("output/step_" + std::to_string( step_num ) + "/JMS_element_" + std::to_string( this->number) + ".txt", lines_);

    }

}