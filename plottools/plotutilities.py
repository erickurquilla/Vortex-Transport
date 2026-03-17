import numpy as np

################################################################################
# Lagrange basis functions
################################################################################

def lagrange_basis_reference_space(p, coords_ref_spa):
    # Coordinates in reference space (xi, eta)
    xi, eta = coords_ref_spa[0], coords_ref_spa[1]
    
    # phi(xi, eta) for each interior node in the element
    phi = [0.0] * ((p + 1) * (p + 2) // 2)

    if p == 0:
        phi[0] = 1.0
    elif p == 1:
        phi[0] = 1.0 - xi - eta
        phi[1] = xi
        phi[2] = eta
    elif p == 2:
        phi[0] = 1.0 - 3.0 * xi - 3.0 * eta + 2.0 * xi * xi + 4.0 * xi * eta + 2.0 * eta * eta
        phi[1] = 4.0 * xi - 4.0 * xi * xi - 4.0 * xi * eta
        phi[2] = -xi + 2.0 * xi * xi
        phi[3] = 4.0 * eta - 4.0 * xi * eta - 4.0 * eta * eta
        phi[4] = 4.0 * xi * eta
        phi[5] = -eta + 2.0 * eta * eta
    elif p == 3:
        phi[0] = 1.0 - 11.0/2.0*xi - 11.0/2.0*eta + 9.0*xi*xi + 18.0*xi*eta + 9.0*eta*eta - 9.0/2.0*xi*xi*xi - 27.0/2.0*xi*xi*eta - 27.0/2.0*xi*eta*eta - 9.0/2.0*eta*eta*eta
        phi[1] = 9.0*xi - 45.0/2.0*xi*xi - 45.0/2.0*xi*eta + 27.0/2.0*xi*xi*xi + 27.0*xi*xi*eta + 27.0/2.0*xi*eta*eta
        phi[2] = -9.0/2.0*xi + 18.0*xi*xi + 9.0/2.0*xi*eta - 27.0/2.0*xi*xi*xi -27.0/2.0*xi*xi*eta
        phi[3] = xi - 9.0/2.0*xi*xi + 9.0/2.0*xi*xi*xi
        phi[4] = 9.0*eta - 45.0/2.0*xi*eta - 45.0/2.0*eta*eta + 27.0/2.0*xi*xi*eta + 27.0*xi*eta*eta + 27.0/2.0*eta*eta*eta
        phi[5] = 27.0*xi*eta - 27.0*xi*xi*eta - 27.0*xi*eta*eta
        phi[6] = -9.0/2.0*xi*eta + 27.0/2.0*xi*xi*eta
        phi[7] = -9.0/2.0*eta + 9.0/2.0*xi*eta + 18.0*eta*eta - 27.0/2.0*xi*eta*eta - 27.0/2.0*eta*eta*eta
        phi[8] = -9.0/2.0*xi*eta + 27.0/2.0*xi*eta*eta
        phi[9] = eta - 9.0/2.0*eta*eta + 9.0/2.0*eta*eta*eta
    else:
        print("ERROR: Unsupported p order for lagrange_basis_reference_space")
        exit(EXIT_FAILURE)

    # Return the value of all Lagrange basis function at position (xi, eta) in reference space
    # First index represents the node number, running between 0 and (p + 1) * (p + 2) / 2
    # The value contained in the Lagrange basis function evaluated at point (xi, eta) in reference space
    return phi

################################################################################
# Gauss quadrature integration
################################################################################

x = [
    0.333333333333333, 0.333333333333333, 0.028844733232685, 0.485577633383657,
    0.485577633383657, 0.485577633383657, 0.485577633383657, 0.028844733232685,
    0.781036849029926, 0.109481575485037, 0.109481575485037, 0.109481575485037,
    0.109481575485037, 0.781036849029926, 0.141707219414880, 0.307939838764121,
    0.307939838764121, 0.550352941820999, 0.550352941820999, 0.141707219414880,
    0.307939838764121, 0.141707219414880, 0.550352941820999, 0.307939838764121,
    0.141707219414880, 0.550352941820999, 0.025003534762686, 0.246672560639903,
    0.246672560639903, 0.728323904597411, 0.728323904597411, 0.025003534762686,
    0.246672560639903, 0.025003534762686, 0.728323904597411, 0.246672560639903,
    0.025003534762686, 0.728323904597411, 0.009540815400299, 0.066803251012200,
    0.066803251012200, 0.923655933587500, 0.923655933587500, 0.009540815400299,
    0.066803251012200, 0.009540815400299, 0.923655933587500, 0.066803251012200,
    0.009540815400299, 0.923655933587500
]

w = [
    0.045408995191377, 0.018362978878233, 0.018362978878233, 0.018362978878233,
    0.022660529717764, 0.022660529717764, 0.022660529717764, 0.036378958422710,
    0.036378958422710, 0.036378958422710, 0.036378958422710, 0.036378958422710,
    0.036378958422710, 0.014163621265528, 0.014163621265528, 0.014163621265528,
    0.014163621265528, 0.014163621265528, 0.014163621265528, 0.004710833481867,
    0.004710833481867, 0.004710833481867, 0.004710833481867, 0.004710833481867,
    0.004710833481867
]

gauss_quadrature = np.zeros( ( len(w) , 2 ) )

for i in range(len(w)):
    gauss_quadrature[i] = [x[2 * i], x[2 * i + 1]]

################################################################################
# READ GRID FILE
################################################################################

def read_input_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines
                key, value = line.strip().split('=')
                if key.startswith('coordinate'):
                    data.append(float(value))
                elif key in ['element_number', 'type', 'right_element', 'left_element', 'vertical_element']:
                    data.append(int(value))
    return data