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


def is_face_inverted(vertices, face):
    A = Matrix(vertices[face[0]])
    B = Matrix(vertices[face[1]])
    C = Matrix(vertices[face[2]])

    edge1 = B - A
    edge2 = C - A

    normal_vector = edge1.cross(edge2)

    centroid = (A + B + C) / 3

    if normal_vector > 0:
        print("face is inverted")
        return False 
    else:
        print("face is not inverted")
        return True 


if __name__ == "__main__":

    x, y, z = symbols('x y z')

    A = Matrix(V[F[0][0]])
    B = Matrix(V[F[0][1]])
    C = Matrix(V[F[0][2]])
    D = Matrix([x, y, z])

    # QUESTION 1
    if args.test == 1:
        #QUESTION 1

        # debugging print statements
        # print("first vector is: ", V[F[0][0]])
        # print("second vector is: ", V[F[0][1]])
        # print("third vector is: ", V[F[0][2]])

        vol_func = lambdify((x, y, z), abs((A.dot(B.cross(C))/6)))

        tetrahedron_vol = vol_func(0,0,0)

        mass = tetrahedron_vol * args.density
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

        vol_func = lambdify((x, y, z), abs((A.dot(B.cross(C))/6)))

        tetrahedron_vol = vol_func(0,0,0)

        rho = args.density

        #Using matrices from question 1
        integrand = (x*A + y*B + z*C + (1-(x + y + z))*D2)*rho

        q2_result = 6*tetrahedron_vol*integrate(integrand, (z, 0, 1-x-y), (y, 0, 1-x), (x, 0, 1))

        mass = tetrahedron_vol * args.density
        print("weighted com = ", q2_result)

        

        
    # QUESTION 3
    elif args.test == 3:
        # debugging print statements
        # print("first vector is: ", V[F[0][0]])
        # print("second vector is: ", V[F[0][1]])
        # print("third vector is: ", V[F[0][2]])

        vol_func = lambdify((x, y, z), abs((A.dot(B.cross(C))/6)))

        tetrahedron_vol = vol_func(0,0,0)

        rho = args.density

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

    x, y, z = symbols('x y z')
    D = Matrix([x, y, z])

    total_vol = 0
    for f in F: 
        vol_func = lambdify((x, y, z), (B-A).dot((C-A).cross(D-A))/6)
        tetrahedron_vol = vol_func(0,0,0)
        if (is_face_inverted(V, f)):  
            total_vol -= tetrahedron_vol
        else:
            total_vol += tetrahedron_vol


    print("Volume of bunny = ", total_vol)

