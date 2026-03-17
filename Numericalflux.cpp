#include <vector>
#include <iostream>
#include <cmath>

#include "Numericalflux.H"

// compute numerical flux given the state vector u in the left and right and the normal vector
std::vector<double> numerical_flux(const std::vector<double>& u_left, const std::vector<double>& u_right, const std::vector<double>& normal_vector){

    // function [F, smag] = flux(UL, UR, n)
    // % PURPOSE: This function calculates the flux for the Euler equations
    // % using the Roe flux function
    // %
    // % INPUTS:
    // %    UL: conservative state vector in left cell
    // %    UR: conservative state vector in right cell
    // %     n: normal pointing from the left cell to the right cell
    // %
    // % OUTPUTS:
    // %  F   : the flux out of the left cell (into the right cell)
    // %  smag: the maximum propagation speed of disturbance
    // %

    double gamma = 1.4;
    double gmi = gamma - 1.0;

    // process left state
    double rL = u_left[0];
    double uL = u_left[1] / rL;
    double vL = u_left[2] / rL;
    double unL = uL * normal_vector[0] + vL * normal_vector[1];
    double qL = pow( pow( u_left[1] , 2 ) + pow( u_left[2] , 2 ) , 0.5 ) / rL;
    double pL = ( gamma - 1 ) * ( u_left[3] - 0.5 * rL * pow( qL , 2 ) );

    if ( ( pL <= 0 ) || ( rL <=0 ) ){
        std::cerr << "1. Error: Non-physical state! " << std::endl;
        std::cout << " pL : "<< pL << " , rL : " << rL << ", not satified ( pL <= 0 ) || ( rL <=0 )" << std::endl;
        std::cout << " u_left[0] = "<<  u_left[0] << std::endl;
        std::cout << " u_left[1] = "<<  u_left[1] << std::endl;
        std::cout << " u_left[2] = "<<  u_left[2] << std::endl;
        std::cout << " u_left[3] = "<<  u_left[3] << std::endl;
        std::cout << " u_right[0] = "<<  u_right[0] << std::endl;
        std::cout << " u_right[1] = "<<  u_right[1] << std::endl;
        std::cout << " u_right[2] = "<<  u_right[2] << std::endl;
        std::cout << " u_right[3] = "<<  u_right[3] << std::endl;
        exit(EXIT_FAILURE);
    }
    
    double rHL = u_left[3] + pL;
    double HL = rHL / rL;
    // double cL = pow( gamma * pL /rL , 0.5 );

    // left flux
    std::vector<double> FL(4);
    FL[0] = rL * unL;
    FL[1] = u_left[1] * unL + pL * normal_vector[0];
    FL[2] = u_left[2] * unL + pL * normal_vector[1];
    FL[3] = rHL * unL;

    // process right state
    double rR = u_right[0];
    double uR = u_right[1] / rR;
    double vR = u_right[2] / rR;
    double unR = uR * normal_vector[0] + vR * normal_vector[1];
    double qR = pow( pow( u_right[1] , 2 ) + pow( u_right[2] , 2 ) , 0.5 ) / rR;
    double pR = ( gamma - 1 ) * ( u_right[3] - 0.5 * rR * pow( qR , 2 ) );

    if ( ( pR <= 0 ) || ( rR <=0 ) ){
        std::cerr << "2. Error: Non-physical state! " << std::endl;
        std::cout << " pR : "<< pR << " , rR : " << rR << ", not satified ( pR <= 0 ) || ( rR <=0 )" << std::endl;
        exit(EXIT_FAILURE);
    }

    double rHR = u_right[3] + pR;
    double HR = rHR / rR;
    // double cR = pow( gamma * pR / rR , 0.5 );

    // right flux
    std::vector<double> FR(4);
    FR[0] = rR * unR;
    FR[1] = u_right[1] * unR + pR * normal_vector[0];
    FR[2] = u_right[2] * unR + pR * normal_vector[1];
    FR[3] = rHR * unR;

    // difference in states
    std::vector<double> du(4);
    du[0] = u_right[0] - u_left[0];
    du[1] = u_right[1] - u_left[1];
    du[2] = u_right[2] - u_left[2];
    du[3] = u_right[3] - u_left[3];

    // Roe average
    double di     = pow( rR / rL , 0.5 );
    double d1     = 1.0 / ( 1.0 + di );

    double ui     = ( di * uR + uL) * d1;
    double vi     = ( di * vR + vL) * d1;
    double Hi     = ( di * HR + HL) * d1;

    double af     = 0.5 * ( ui * ui + vi * vi );
    double ucp    = ui * normal_vector[0] + vi * normal_vector[1];
    double c2     = gmi * ( Hi - af );
    
    if ( c2 <= 0 ){
        std::cerr << "3. Error: Non-physical state! " << std::endl;
        std::cout << " c2 : "<< pR << " , not satified ( c2 <= 0 )" << std::endl;
        exit(EXIT_FAILURE);
    }

    double ci     = pow( c2 , 2 );
    double ci1    = 1.0 / ci;

    // eigenvalues
    // z = zeros(3,1);
    std::vector<double> l = {0,0,0};
    l[0] = ucp + ci ;
    l[1] = ucp - ci ;
    l[2] = ucp;

    // entropy fix
    double epsilon = ci * 0.1;
    for (int i = 0; i < 3; ++i) {
        if ( ( l[i] < epsilon ) && ( l[i] > -1 * epsilon ) ){
            l[i] = 0.5 * (epsilon + l[i] * l[i] / epsilon);
        }
    }

    for (int i = 0; i < 3; ++i) {
        l[i] = std::abs(l[i]);
    }
    double l3 = l[2];

    // Average and half-difference of 1st and 2nd eigenvalues
    double s1 = 0.5 * (l[0] + l[1]);
    double s2 = 0.5 * (l[0] - l[1]);

    // Left eigenvector product generators
    double G1 = gmi * (af * du[0] - ui * du[1] - vi * du[2] + du[3]);
    double G2 = -1 * ucp * du[0] + du[1] * normal_vector[0] + du[2] * normal_vector[1];

    // Required functions of G1 and G2
    double C1 = G1 * (s1 - l3) * ci1 * ci1 + G2 * s2 * ci1;
    double C2 = G1 * s2 * ci1 + G2 * (s1 - l3);

    // Flux assembly
    std::vector<double> F(4);
    F[0] = 0.5 * (FL[0] + FR[0]) - 0.5 * (l3 * du[0] + C1);
    F[1] = 0.5 * (FL[1] + FR[1]) - 0.5 * (l3 * du[1] + C1 * ui + C2 * normal_vector[0]);
    F[2] = 0.5 * (FL[2] + FR[2]) - 0.5 * (l3 * du[2] + C1 * vi + C2 * normal_vector[1]);
    F[3] = 0.5 * (FL[3] + FR[3]) - 0.5 * (l3 * du[3] + C1 * Hi + C2 * ucp);

    return F;
}