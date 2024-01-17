from sympy import *
import numpy as np
import igl
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--file", type=str, default="cube.obj")
parser.add_argument("--density", type=float, default=1000, help="density in kg/m^3")
parser.add_argument("--scale", nargs=3, type=float, default=(1,1,1),help="scale the mesh")
parser.add_argument("--translate", nargs=3, type=float, default=(0,0,0),help="translate (after scale)")
parser.add_argument("--test", type=int, help="run a numbered unit test")
args = parser.parse_args()
V, _, _, F, _, _ = igl.read_obj(args.file)
V = V * args.scale
V = V + args.translate


if __name__ == "__main__":

    x, y, z = symbols('x y z')

    A = Matrix(V[F[0][0]])
    B = Matrix(V[F[0][1]])
    C = Matrix(V[F[0][2]])
    D = Matrix([x, y, z])

    vol_func = lambdify((x, y, z), abs((A.dot(B.cross(C))/6)))

    tetrahedron_vol = vol_func(0,0,0)

    rho = args.density

    # QUESTION 1
    if args.test == 1:
        #QUESTION 1

        # debugging print statements
        # print("first vector is: ", V[F[0][0]])
        # print("second vector is: ", V[F[0][1]])
        # print("third vector is: ", V[F[0][2]])

        mass = tetrahedron_vol * rho
        print("vol = ", tetrahedron_vol)
        print("mass = ", mass)



    # QUESTION 2
    elif args.test == 2:
        #QUESTION 1
        D2 = Matrix([0, 0, 0])

        # debugging print statements
        # print("first vector is: ", V[F[0][0]])
        # print("second vector is: ", V[F[0][1]])
        # print("third vector is: ", V[F[0][2]])

        rho = args.density

        #Using matrices from question 1
        integrand = (x*A + y*B + z*C + (1-(x + y + z))*D2)*rho

        q2_result = 6*tetrahedron_vol*integrate(integrand, (z, 0, 1-x-y), (y, 0, 1-x), (x, 0, 1))

        print("weighted com = ", q2_result)

        

        

        
    # QUESTION 3
    elif args.test == 3:
        # debugging print statements
        # print("first vector is: ", V[F[0][0]])
        # print("second vector is: ", V[F[0][1]])
        # print("third vector is: ", V[F[0][2]])

        p = x*A + y*B + z*C

        r_matrix = Matrix([[0, -p[2], p[1]],[p[2], 0, -p[0]], [-p[1], p[0], 0]])

        #debugging tools

        #pprint(r_matrix)

        #pprint(r_matrix.transpose())

        inertia_integrand = r_matrix.transpose()*r_matrix

        j = 6*tetrahedron_vol*rho*integrate(inertia_integrand, (z, 0, 1 - x - y), (y, 0, 1 - x), (x, 0, 1))

        print("J = ")
        pprint(j)

    #QUESTION 4

    A_1, B_1, C_1, A_2, B_2, C_2, A_3, B_3, C_3 = symbols('A_1 B_1 C_1 A_2 B_2 C_2 A_3 B_3 C_3')
    A = Matrix([A_1, A_2, A_3])
    B = Matrix([B_1, B_2, B_3])
    C = Matrix([C_1, C_2, C_3])
    D = Matrix([x, y, z])
    
    total_vol_q4 = 0
    com_q4 = Matrix([0, 0, 0])
    j_q4 = Matrix([[0, 0, 0],[0, 0, 0], [0, 0, 0]])

    vol_tetra_q4 = A.dot(B.cross(C))/6
    integrand_com_q4 = (x*A + y*B + z*C)*rho

    p_q4 = x*A + y*B + z*C
    r_matrix_q4 = Matrix([[0, -p_q4[2], p_q4[1]],[p_q4[2], 0, -p_q4[0]], [-p_q4[1], p_q4[0], 0]])
    inertia_integrand_q4 = r_matrix_q4.transpose()*r_matrix_q4

    vol_func_q4 = lambdify(([A_1, A_2, A_3], [B_1, B_2, B_3], [C_1, C_2, C_3]), vol_tetra_q4, 'sympy')
    weighted_com_func_q4 = lambdify(([A_1, A_2, A_3], [B_1, B_2, B_3], [C_1, C_2, C_3]), integrate(integrand_com_q4, (z, 0, 1-x-y), (y, 0, 1-x), (x, 0, 1)), 'sympy')
    j_func_q4 = lambdify(([A_1, A_2, A_3], [B_1, B_2, B_3], [C_1, C_2, C_3]), integrate(inertia_integrand_q4, (z, 0, 1 - x - y), (y, 0, 1 - x), (x, 0, 1)), 'sympy')

    j_q4 = Matrix([[0, 0, 0],[0, 0, 0], [0, 0, 0]])
    for f in F: 
        A = Matrix(V[f[0]])
        B = Matrix(V[f[1]])
        C = Matrix(V[f[2]])

        tetrahedron_vol_q4 = vol_func_q4(V[f[0]],V[f[1]],V[f[2]])
        total_vol_q4 += tetrahedron_vol_q4

        #weighted com
        com_q4 += weighted_com_func_q4(V[f[0]],V[f[1]],V[f[2]])*tetrahedron_vol_q4*6

        #rotational intertia
        j_q4 += j_func_q4(V[f[0]],V[f[1]],V[f[2]])*tetrahedron_vol_q4*6*rho

        # p_q4 = x*A_q4 + y*B_q4 + z*C_q4
        # r_matrix_q4 = Matrix([[0, -p_q4[2], p_q4[1]],[p_q4[2], 0, -p_q4[0]], [-p_q4[1], p_q4[0], 0]])
        # inertia_integrand_q4 = r_matrix_q4.transpose()*r_matrix_q4
        # j_q4 += 6*tetrahedron_vol_q4*rho*integrate(inertia_integrand_q4, (z, 0, 1 - x - y), (y, 0, 1 - x), (x, 0, 1))

    mass_q4 = total_vol_q4 * args.density

    print("volume = ", np.round(total_vol_q4, decimals = 3))
    print("mass = ", np.round(mass_q4, decimals = 3))
    print("com = ", com_q4/mass_q4)

    numpy_ar = np.array(j_q4, dtype=float)
    rounded_j_q4_ar = np.round(numpy_ar, decimals = 3)
    j_q4_rounded = Matrix(rounded_j_q4_ar)
    print("J = ")
    pprint(j_q4_rounded)


