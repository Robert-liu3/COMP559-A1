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


#QUESTION 1

x, y, z = symbols('x y z')

A = Matrix(V[F[0][0]])
B = Matrix(V[F[0][1]])
C = Matrix(V[F[0][2]])
D = Matrix([x, y, z])

# debugging print statements
# print("first vector is: ", V[F[0][0]])
# print("second vector is: ", V[F[0][1]])
# print("third vector is: ", V[F[0][2]])

vol_func = lambdify((x, y, z), abs(C.dot(A.cross(B))/6), "numpy")

tetrahedron_vol = vol_func(0,0,0)

mass = tetrahedron_vol * args.density


#QUESTION 2

D2 = Matrix([0, 0, 0])

rho = args.density

#Using matrices from question 1
integrand = (x*A + y*B + z*C + (1-(x + y + z))*D2)*rho

q2_result = 6*tetrahedron_vol*integrate(integrand, (z, 0, 1-x-y), (y, 0, 1-x), (x, 0, 1))


#QUESTION 3

# r_A = A - D2

p = x*A + y*B + z*C + (1-(x + y + z))*D2

r_matrix = Matrix([[0, -p[2], p[1]],[p[2], 0, -p[0]], [-p[1], p[0], 0]])

#debugging tools

pprint(r_matrix)

pprint(r_matrix.transpose())

inertia_integrand = r_matrix.transpose()*r_matrix

j = rho*integrate(inertia_integrand, (z, 0, 1 - x - y), (y, 0, 1 - x), (x, 0, 1))







# QUESTION 1
if args.test == 1:
    print("vol = ", tetrahedron_vol)


    print("mass = ", mass)



# QUESTION 2
elif args.test == 2:

    print("weighted com = ", q2_result)

    
# QUESTION 3
elif args.test == 3:

    print("J = ", j)
    pprint(j)

