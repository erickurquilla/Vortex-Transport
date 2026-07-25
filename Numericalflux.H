#ifndef NUMERICALFLUX_H
#define NUMERICALFLUX_H

// compute numerical flux given the state vector u in the left and right and the normal vector
// the result is written into `flux`, which must have size 4 (no allocation happens here:
// this function sits in the innermost hot loop, one call per quadrature point per side)
void numerical_flux(const std::vector<double>& u_left, const std::vector<double>& u_right, const std::vector<double>& normal_vector, std::vector<double>& flux);

#endif // NUMERICALFLUX_H
