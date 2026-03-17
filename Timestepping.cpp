#include <vector>
#include <cmath>
#include <iostream>

#include "Element.H"
#include "Evolve.H"

// Compute dynamical quantites of the DG method
void evolve_elem(Evolve_element* evo_elemts, int n_elements) {
    // Compute require quantities for time evolution
    for (int n = 0; n < n_elements ; ++n) {
        evo_elemts[n].compute_U_plus_minus();      // compute U on the element boundaries
        evo_elemts[n].compute_numerical_flux();    // compute numerical flux
        evo_elemts[n].integrate_numerical_flux();  // compute the vector result of integrating the numerical flux
        evo_elemts[n].compute_stiffness_vector();  // compute sitffness vector
        evo_elemts[n].compute_residual_vector();   // compute residual vector
        evo_elemts[n].compute_time_derivative_U(); // compute time derivative of U
    }
}

// time stepping with forward euler
void forward_euler(Element* elemts, Evolve_element* evo_elemts, double& t_step, const int& n_elements){

    // Compute require quantities for time evolution
    evolve_elem(evo_elemts, n_elements);

    // define hidrodynamics quantities
    double rho, u, v, E, H, p;
    double gamma = 1.4;

    // loop over elements
    for (int n = 0; n < n_elements ; ++n) {
        // loop over all the interior nodes of this element
        for (int i = 0; i < ( elemts[n].p + 1 ) * ( elemts[n].p + 2 ) / 2; ++i) {
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                elemts[n].hidrodynamics_vector_U[i][k] += t_step * evo_elemts[n].DG_time_derivative_U[i][k]; 
            }

            rho = elemts[n].hidrodynamics_vector_U[i][0]; // density
            u   = elemts[n].hidrodynamics_vector_U[i][1] / rho; // horizontal velocity
            v   = elemts[n].hidrodynamics_vector_U[i][2] / rho; // vertical velocity
            E   = elemts[n].hidrodynamics_vector_U[i][3] / rho; // energy
            p   = rho * ( gamma - 1 ) * ( E - ( pow( u , 2) + pow( v , 2) ) / 2 ); // pressure
            H   = E + p / rho; // Entalpy

            // initialize the hidrodinamic vector f, x component  
            elemts[n].hidrodynamics_vector_F[i][0][0] = rho * u;
            elemts[n].hidrodynamics_vector_F[i][0][1] = rho * pow( u , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][0][2] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][0][3] = rho * u * H;

            // initialize the hidrodinamic vector f, y component  
            elemts[n].hidrodynamics_vector_F[i][1][0] = rho * v;
            elemts[n].hidrodynamics_vector_F[i][1][1] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][1][2] = rho * pow( v , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][1][3] = rho * v * H;
        }
        
        // advance in time
        elemts[n].time += t_step;
    }
}

// time stepping with rk4
void rk4(Element* elemts, Evolve_element* evo_elemts, double& t_step, const int& n_elements){

    // RK4 vector k1, k1, k3 and k4
    // The first index runs over elements
    // The second index runs over interior nodes
    // The third inxed runs between 0 and 4 for hidrodynamics quantities: 0: rho, 1: rho u, 2: rho v, 3: rho E.
    std::vector<std::vector<std::vector<double>>> k1(n_elements, std::vector<std::vector<double>>((elemts[0].p + 1) * (elemts[0].p + 2) / 2, std::vector<double>(4)));
    std::vector<std::vector<std::vector<double>>> k2(n_elements, std::vector<std::vector<double>>((elemts[0].p + 1) * (elemts[0].p + 2) / 2, std::vector<double>(4)));
    std::vector<std::vector<std::vector<double>>> k3(n_elements, std::vector<std::vector<double>>((elemts[0].p + 1) * (elemts[0].p + 2) / 2, std::vector<double>(4)));
    std::vector<std::vector<std::vector<double>>> k4(n_elements, std::vector<std::vector<double>>((elemts[0].p + 1) * (elemts[0].p + 2) / 2, std::vector<double>(4)));

    // Hidrodynamic vector U
    // The first index runs over elements
    // The second index runs over interior nodes
    // The third inxed runs between 0 and 4 for hidrodynamics quantities: 0: rho, 1: rho u, 2: rho v, 3: rho E.
    std::vector<std::vector<std::vector<double>>> U(n_elements, std::vector<std::vector<double>>((elemts[0].p + 1) * (elemts[0].p + 2) / 2, std::vector<double>(4)));

    // loop over elements
    for (int n = 0; n < n_elements ; ++n) {
        // loop over all the interior nodes of this element
        for (int i = 0; i < ( elemts[n].p + 1 ) * ( elemts[n].p + 2 ) / 2; ++i) {
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                // create a copy of hydrodynamic vector U                
                U[n][i][k] = elemts[n].hidrodynamics_vector_U[i][k];
            }
        }
    }

    // define useful hidrodynamics quantities
    double rho, u, v, E, H, p;
    double gamma = 1.4;

    //
    // COMPUTE K1
    //

    // Compute require quantities for time evolution
    evolve_elem(evo_elemts, n_elements);

    // loop over elements
    for (int n = 0; n < n_elements ; ++n) {
        // loop over all the interior nodes of this element
        for (int i = 0; i < ( elemts[n].p + 1 ) * ( elemts[n].p + 2 ) / 2; ++i) {
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                
                // k1 = h * f( y_n )
                k1[n][i][k] = t_step * evo_elemts[n].DG_time_derivative_U[i][k];
                // y_n -> y_n + k1/2
                elemts[n].hidrodynamics_vector_U[i][k] = U[n][i][k] + k1[n][i][k] / 2.0; 

            }

            rho = elemts[n].hidrodynamics_vector_U[i][0]; // density
            u   = elemts[n].hidrodynamics_vector_U[i][1] / rho; // horizontal velocity
            v   = elemts[n].hidrodynamics_vector_U[i][2] / rho; // vertical velocity
            E   = elemts[n].hidrodynamics_vector_U[i][3] / rho; // energy
            p   = rho * ( gamma - 1 ) * ( E - ( pow( u , 2) + pow( v , 2) ) / 2 ); // pressure
            H   = E + p / rho; // Entalpy

            // initialize the hidrodinamic vector f, x component  
            elemts[n].hidrodynamics_vector_F[i][0][0] = rho * u;
            elemts[n].hidrodynamics_vector_F[i][0][1] = rho * pow( u , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][0][2] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][0][3] = rho * u * H;

            // initialize the hidrodinamic vector f, y component  
            elemts[n].hidrodynamics_vector_F[i][1][0] = rho * v;
            elemts[n].hidrodynamics_vector_F[i][1][1] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][1][2] = rho * pow( v , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][1][3] = rho * v * H;
        }
    }

    //
    // COMPUTE K2
    //

    // Compute require quantities for time evolution
    evolve_elem(evo_elemts, n_elements);

    // loop over elements
    for (int n = 0; n < n_elements ; ++n) {
        // loop over all the interior nodes of this element
        for (int i = 0; i < ( elemts[n].p + 1 ) * ( elemts[n].p + 2 ) / 2; ++i) {
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                
                // k2 = h * f( y_n + k1/2 )
                k2[n][i][k] = t_step * evo_elemts[n].DG_time_derivative_U[i][k];
                // y_n -> y_n + k2/2
                elemts[n].hidrodynamics_vector_U[i][k] = U[n][i][k] + k2[n][i][k] / 2.0; 

            }

            rho = elemts[n].hidrodynamics_vector_U[i][0]; // density
            u   = elemts[n].hidrodynamics_vector_U[i][1] / rho; // horizontal velocity
            v   = elemts[n].hidrodynamics_vector_U[i][2] / rho; // vertical velocity
            E   = elemts[n].hidrodynamics_vector_U[i][3] / rho; // energy
            p   = rho * ( gamma - 1 ) * ( E - ( pow( u , 2) + pow( v , 2) ) / 2 ); // pressure
            H   = E + p / rho; // Entalpy

            // initialize the hidrodinamic vector f, x component  
            elemts[n].hidrodynamics_vector_F[i][0][0] = rho * u;
            elemts[n].hidrodynamics_vector_F[i][0][1] = rho * pow( u , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][0][2] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][0][3] = rho * u * H;

            // initialize the hidrodinamic vector f, y component  
            elemts[n].hidrodynamics_vector_F[i][1][0] = rho * v;
            elemts[n].hidrodynamics_vector_F[i][1][1] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][1][2] = rho * pow( v , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][1][3] = rho * v * H;
        }
    }

    //
    // COMPUTE K3
    //

    // Compute require quantities for time evolution
    evolve_elem(evo_elemts, n_elements);

    // loop over elements
    for (int n = 0; n < n_elements ; ++n) {
        // loop over all the interior nodes of this element
        for (int i = 0; i < ( elemts[n].p + 1 ) * ( elemts[n].p + 2 ) / 2; ++i) {
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                
                // k3 = h * f( y_n + k2/2 )
                k3[n][i][k] = t_step * evo_elemts[n].DG_time_derivative_U[i][k];
                // y_n -> y_n + k3          
                elemts[n].hidrodynamics_vector_U[i][k] = U[n][i][k] + k3[n][i][k];
                
            }

            rho = elemts[n].hidrodynamics_vector_U[i][0]; // density
            u   = elemts[n].hidrodynamics_vector_U[i][1] / rho; // horizontal velocity
            v   = elemts[n].hidrodynamics_vector_U[i][2] / rho; // vertical velocity
            E   = elemts[n].hidrodynamics_vector_U[i][3] / rho; // energy
            p   = rho * ( gamma - 1 ) * ( E - ( pow( u , 2) + pow( v , 2) ) / 2 ); // pressure
            H   = E + p / rho; // Entalpy

            // initialize the hidrodinamic vector f, x component  
            elemts[n].hidrodynamics_vector_F[i][0][0] = rho * u;
            elemts[n].hidrodynamics_vector_F[i][0][1] = rho * pow( u , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][0][2] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][0][3] = rho * u * H;

            // initialize the hidrodinamic vector f, y component  
            elemts[n].hidrodynamics_vector_F[i][1][0] = rho * v;
            elemts[n].hidrodynamics_vector_F[i][1][1] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][1][2] = rho * pow( v , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][1][3] = rho * v * H;
        }
    }

    //
    // COMPUTE K4
    //

    // Compute require quantities for time evolution
    evolve_elem(evo_elemts, n_elements);

    // loop over elements
    for (int n = 0; n < n_elements ; ++n) {
        // loop over all the interior nodes of this element
        for (int i = 0; i < ( elemts[n].p + 1 ) * ( elemts[n].p + 2 ) / 2; ++i) {
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                
                // k4 = h * f( y_n + k3 )
                k4[n][i][k] = t_step * evo_elemts[n].DG_time_derivative_U[i][k];

            }
        }
    }

    //
    // COMPUTE NEW U AND F
    //

    // loop over elements
    for (int n = 0; n < n_elements ; ++n) {
        // loop over all the interior nodes of this element
        for (int i = 0; i < ( elemts[n].p + 1 ) * ( elemts[n].p + 2 ) / 2; ++i) {
            // loop over hidrodynamics indices
            for (int k = 0; k < 4; ++k) {
                
                // y_n+1 = y_n + 1/6 ( k1 + 2 * k2 + 2 * k3 + k4 )          
                elemts[n].hidrodynamics_vector_U[i][k] = U[n][i][k] + ( 1.0 / 6.0 ) * ( k1[n][i][k] + 2 * k2[n][i][k] + 2 * k3[n][i][k] + k4[n][i][k] ); 

            }

            rho = elemts[n].hidrodynamics_vector_U[i][0]; // density
            u   = elemts[n].hidrodynamics_vector_U[i][1] / rho; // horizontal velocity
            v   = elemts[n].hidrodynamics_vector_U[i][2] / rho; // vertical velocity
            E   = elemts[n].hidrodynamics_vector_U[i][3] / rho; // energy
            p   = rho * ( gamma - 1 ) * ( E - ( pow( u , 2) + pow( v , 2) ) / 2 ); // pressure
            H   = E + p / rho; // Entalpy

            // initialize the hidrodinamic vector f, x component  
            elemts[n].hidrodynamics_vector_F[i][0][0] = rho * u;
            elemts[n].hidrodynamics_vector_F[i][0][1] = rho * pow( u , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][0][2] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][0][3] = rho * u * H;

            // initialize the hidrodinamic vector f, y component  
            elemts[n].hidrodynamics_vector_F[i][1][0] = rho * v;
            elemts[n].hidrodynamics_vector_F[i][1][1] = rho * u * v;
            elemts[n].hidrodynamics_vector_F[i][1][2] = rho * pow( v , 2 ) + p;
            elemts[n].hidrodynamics_vector_F[i][1][3] = rho * v * H;
        }

        // advance in time
        elemts[n].time += t_step; 
    }
}