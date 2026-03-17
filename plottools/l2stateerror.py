import pandas as pd
import glob
import numpy as np

from plotutilities import lagrange_basis_reference_space, gauss_quadrature, w, read_input_file

################################################################################
# L2 STATE ERROR
################################################################################

# Get directory names
dir_names = glob.glob("output/step_*")
dir_names = sorted(dir_names, key=lambda x: int(x.lstrip("output/step_")))

# t=0.0 : First data file generated
file_names_i = glob.glob(dir_names[0] + "/element_*.txt")
file_names_i = sorted(file_names_i, key=lambda x: int(x.split("/element_")[1].rstrip(".txt")))

# t=14.1421 : Last data file generated
file_names_f = glob.glob(dir_names[-1] + "/element_*.txt")
file_names_f = sorted(file_names_f, key=lambda x: int(x.split("/element_")[1].rstrip(".txt")))

# Get grid file names matching the pattern "element*.txt"
grid_file_names = glob.glob("grid/element_*.txt")
grid_file_names = sorted(grid_file_names, key=lambda x: int(x.lstrip("grid/element_").rstrip(".txt")))

# save integral over space of hidrodynamic vector u
integral = np.zeros(4)

# runs over elements
for i in range(len(file_names_i)):

    # read grid data
    data = read_input_file(grid_file_names[i])
    # save elements vertex
    vertices_coords_phys_space =  [ [ data[1] , data[2] ] , [ data[3] , data[4] ] , [ data[5] , data[6] ] ]
    # compute jacobian
    jacobian = [
        [vertices_coords_phys_space[1][0] - vertices_coords_phys_space[0][0], vertices_coords_phys_space[2][0] - vertices_coords_phys_space[0][0]],
        [vertices_coords_phys_space[1][1] - vertices_coords_phys_space[0][1], vertices_coords_phys_space[2][1] - vertices_coords_phys_space[0][1]]
    ]
    # Compute determinant of jacobian
    determinant_jacobian = jacobian[0][0] * jacobian[1][1] - jacobian[1][0] * jacobian[0][1]

    # read data at t=0.0
    df_i = pd.read_csv(file_names_i[i], sep='\s+')

    x_ele_i = np.array(df_i['x'])
    y_ele_i = np.array(df_i['y'])
    u0_ele_i = np.array(df_i['u0'])
    u1_ele_i = np.array(df_i['u1'])
    u2_ele_i = np.array(df_i['u2'])
    u3_ele_i = np.array(df_i['u3'])

    # read data at t=14.1421
    df_f = pd.read_csv(file_names_f[i], sep='\s+')

    x_ele_f = np.array(df_f['x'])
    y_ele_f = np.array(df_f['y'])
    u0_ele_f = np.array(df_f['u0'])
    u1_ele_f = np.array(df_f['u1'])
    u2_ele_f = np.array(df_f['u2'])
    u3_ele_f = np.array(df_f['u3'])

    # ( delta u ) ^ 2
    del_u0 = (u0_ele_f - u0_ele_i)**2
    del_u1 = (u1_ele_f - u1_ele_i)**2
    del_u2 = (u2_ele_f - u2_ele_i)**2
    del_u3 = (u3_ele_f - u3_ele_i)**2

    p = 3

    basis_in_quad_points = np.zeros( ( int( ( p + 1 ) * ( p + 2) / 2 ) , len( w ) ) )

    # computing basis function in the gauss quadrature points
    for i in range(len(w)):
        basis = lagrange_basis_reference_space( p , gauss_quadrature[i] )
        for j in range( int( ( p + 1 ) * ( p + 2) / 2 ) ):
            basis_in_quad_points[j][i] = basis[j]

    delta_U_in_qp = np.zeros( ( len(w) , 4 ) )

    # intepolating delta u to quadrature points
    for i in range(len(w)):
        for j in range( int( ( p + 1 ) * ( p + 2) / 2 ) ):
            delta_U_in_qp[i][0] += del_u0[j] * basis_in_quad_points[j][i]
            delta_U_in_qp[i][1] += del_u1[j] * basis_in_quad_points[j][i]
            delta_U_in_qp[i][2] += del_u2[j] * basis_in_quad_points[j][i]
            delta_U_in_qp[i][3] += del_u3[j] * basis_in_quad_points[j][i]

    # integrating 
    integral[0] += determinant_jacobian * np.sum( np.array(w) * delta_U_in_qp[:, 0])
    integral[1] += determinant_jacobian * np.sum( np.array(w) * delta_U_in_qp[:, 1])
    integral[2] += determinant_jacobian * np.sum( np.array(w) * delta_U_in_qp[:, 2])
    integral[3] += determinant_jacobian * np.sum( np.array(w) * delta_U_in_qp[:, 3])

# stimating error
error = np.sqrt( ( 1.0 / 25.0 ) * integral )
print(f' L2 state error = {error}')

