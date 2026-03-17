#ifndef TIMESTEPPING_H
#define TIMESTEPPING_H

#include <vector>

#include "Element.H"
#include "Evolve.H"

void forward_euler(Element* elemts, Evolve_element* evo_elemts, double& t_step, const int& n_elements); // time stepping with forward euler
void rk4(Element* elemts, Evolve_element* evo_elemts, double& t_step, const int& n_elements); // time stepping with rk4

#endif // TIMESTEPPING_H