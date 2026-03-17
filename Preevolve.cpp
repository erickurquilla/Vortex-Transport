#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "Lagrangebasis.H"

// compute the inverse of the mass matrix in reference space 
// mass_ij = int in T phi_i phi_j dT
// T is an triangle in reference space with vertex (0,0), (1,0) and (0,1) in reference space.
// return a matrix of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
std::vector<std::vector<double>> inverse_mass_matrix_reference_space(const int& p, const std::vector<std::vector<double>>& gauss_area_int){

    // get size, that is the number of quadrature point for integration
    int size = gauss_area_int.size();

    // this vector store the values of the lagrange polinomial in evaluated in the quadrature points
    // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
    // second item runs over the evaluation of the lagrange poliniam in the quadrature points
    std::vector<std::vector<double>> phi_in_quadrature_points( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( size ) );

    // evaluate the lagrange polinomial in the quadrature points
    for (int i = 0; i < size; ++i) {

        std::vector<double> xi_eta_gauss(2);
        xi_eta_gauss[0] = gauss_area_int[i][0];
        xi_eta_gauss[1] = gauss_area_int[i][1];

        std::vector<double> phi_in_xi_eta_gauss = lagrange_basis_reference_space( p , xi_eta_gauss ); 

        // evaluate the lagrange polinomial in the quadrature points
        for (int j = 0; j < ( p + 1 ) * ( p + 2 ) / 2; ++j) {
            phi_in_quadrature_points[j][i] = phi_in_xi_eta_gauss[j];
        }
    }

    // define mass matrix
    std::vector<std::vector<double>> mass_matrix( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( ( p + 1 ) * ( p + 2 ) / 2 ) );

    // integrate mass matrix using gauss methods
    // mass_ij = int in T phi_i phi_j dT , T is an triangle with vertex (0,0), (1,0) and (0,1) in reference space.
    for (int i = 0; i < ( p + 1 ) * ( p + 2 ) / 2; ++i) {
        for (int j = 0; j < ( p + 1 ) * ( p + 2 ) / 2; ++j) {
            mass_matrix[i][j] = 0;
            for (int n = 0; n < size; ++n) {
                mass_matrix[i][j] += phi_in_quadrature_points[i][n] * phi_in_quadrature_points[j][n] * gauss_area_int[n][2];
            }
        }
    } 

    // Convert the mass matrix to an Eigen matrix just to compute the inverse matrix
    Eigen::MatrixXd eigen_mass_matrix( ( p + 1 ) * ( p + 2 ) / 2 ,  ( p + 1 ) * ( p + 2 ) / 2 );
    for (int i = 0; i <  ( p + 1 ) * ( p + 2 ) / 2 ; ++i) {
        for (int j = 0; j <  ( p + 1 ) * ( p + 2 ) / 2 ; ++j) {
            eigen_mass_matrix(i, j) = mass_matrix[i][j];
        }
    }

    // Compute the inverse of the mass matrix using Eigen
    Eigen::MatrixXd inv_matrix = eigen_mass_matrix.inverse();

    // Convert the inv_matrix to an vector matrix to return the right function type
    std::vector<std::vector<double>> mass_matrix_inverse( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( ( p + 1 ) * ( p + 2 ) / 2 ) );
    for (int i = 0; i <  ( p + 1 ) * ( p + 2 ) / 2 ; ++i) {
        for (int j = 0; j <  ( p + 1 ) * ( p + 2 ) / 2 ; ++j) {
            mass_matrix_inverse[i][j] = inv_matrix(i, j);
        }
    }

    // return inverse mass matrix
    // size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
    return mass_matrix_inverse;
}

// compute stiffness matrix in reference space 
// S_ij = integral in T of ( Nabla phi_i ) phi_j dT
// T is an triangle in reference space with vertex (0,0), (1,0) and (0,1) in reference space.
// return a two dimenional vector in reference space xi and eta, the components of xi and eta are matrices of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2
// first index run between spacial components in reference space. 0: xi and 1 eta.
// second and third index run over matrix inidices of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ]
std::vector<std::vector<std::vector<double>>> sitffness_matrix_reference_space(const int& p, const std::vector<std::vector<double>>& gauss_area_int){

    // get size, that is the number of quadrature point for integration
    int size = gauss_area_int.size();

    // this vector store the values of the lagrange polinomial evaluated in the quadrature points
    // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
    // second index runs quadrature points
    std::vector<std::vector<double>> phi_in_quadrature_points( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( size ) );

    // this vector store the values of the gradients of the lagrange polinomial evaluated in the quadrature points
    // first index runs over interior nodes number, that is, the lagrange polinimial that is one on this node
    // second index runs quadrature points
    // third index runs over the x:0 and y:1 component of the gradient of the lagrange poliniam evaluated at the quadrature points
    std::vector<std::vector<std::vector<double>>> gradient_phi_in_quadrature_points( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<std::vector<double>>( size, std::vector<double>(2)) );

    // evaluate the lagrange polinomial in the quadrature points
    // evaluate the gradient of the lagrange polinomials in the quadrature points
    for (int i = 0; i < size; ++i) {

        std::vector<double> xi_eta_gauss(2);
        xi_eta_gauss[0] = gauss_area_int[i][0];
        xi_eta_gauss[1] = gauss_area_int[i][1];

        std::vector<double> phi_in_xi_eta_gauss = lagrange_basis_reference_space( p , xi_eta_gauss ); 
        std::vector<std::vector<double>> gradient_phi_in_xi_eta_gauss = lagrange_basis_gradient_reference_space( p , xi_eta_gauss ); 

        // save the lagrange polinomials the their gradient in the quadrature points in class variables
        for (int j = 0; j < ( p + 1 ) * ( p + 2 ) / 2; ++j) {
            phi_in_quadrature_points[j][i] = phi_in_xi_eta_gauss[j];
            gradient_phi_in_quadrature_points[j][i][0] = gradient_phi_in_xi_eta_gauss[j][0];
            gradient_phi_in_quadrature_points[j][i][1] = gradient_phi_in_xi_eta_gauss[j][1];
        }
    }

    // define sitffness matrix
    std::vector<std::vector<std::vector<double>>> stiffness_matrix( 2, std::vector<std::vector<double>>( ( p + 1 ) * ( p + 2 ) / 2 , std::vector<double>( ( p + 1 ) * ( p + 2 ) / 2) ) );

    // integrate stiffness matrix using gauss methods
    // S_ij = integral in T of ( Nabla phi_i ) phi_j dT
    for (int i = 0; i < ( p + 1 ) * ( p + 2 ) / 2; ++i) {
        for (int j = 0; j < ( p + 1 ) * ( p + 2 ) / 2; ++j) {
            stiffness_matrix[0][i][j] = 0;
            stiffness_matrix[1][i][j] = 0;
            for (int n = 0; n < size; ++n) {
                stiffness_matrix[0][i][j] += gradient_phi_in_quadrature_points[i][n][0] * phi_in_quadrature_points[j][n] * gauss_area_int[n][2];
                stiffness_matrix[1][i][j] += gradient_phi_in_quadrature_points[i][n][1] * phi_in_quadrature_points[j][n] * gauss_area_int[n][2];
            }
        }
    }

    // return stiffness matrix
    // form :  hat{e}_xi * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ] + hat{e}_eta * matrix[ ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ]     
    // first index run between spacial components in reference space. 0: xi and 1 eta.
    // second and third index run over matrix inidices of size ( p + 1 ) * ( p + 2 ) / 2 by ( p + 1 ) * ( p + 2 ) / 2 ]
    return stiffness_matrix;
}