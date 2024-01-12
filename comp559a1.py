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


# QUESTION 1
x, y, z = symbols('x3 y3 z3')

A = Matrix(V[F[0][0]])
B = Matrix(V[F[0][1]])
C = Matrix(V[F[0][2]])
D = Matrix([x, y, z])

# debugging print statements
# print("first vector is: ", V[F[0][0]])
# print("second vector is: ", V[F[0][1]])
# print("third vector is: ", V[F[0][2]])

vol_func = lambdify((x, y, z), abs(C.dot(A.cross(B))/6), "numpy")

tetrahedron_volume = vol_func(0,0,0)

mass = tetrahedron_volume * args.density



# QUESTION 2



if args.test == 1:

    print("Volume of tetrahedron: ", tetrahedron_volume)


    print("Mass of tetrahedron: ", mass)

elif args.test == 2:

    rho = symbols('rho', positive=True)

    baycentric_func = lambdify((x, y, z), (C.dot(A.cross(B))/6), "numpy")

    


