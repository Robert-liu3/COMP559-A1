from sympy import *
import numpy as np
import igl
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--file", type=str, default="cube.obj")
parser.add_argument("--density", type-float, default=1000, help="density in kg/m^3")
parser.add_argument("--scale", nargs=3, type=float, default=(1,1,1),help="scale the mesh")
parser.add_argument("--translate", nargs=3, type=float, default=(0,0,0),help="translate (after scale)")
parser.add_argument("--test", type=int, help="run a numbered unit test")
args = parser.parse_args()
V, _, _, F, _, _ = igl.read_obj(args.file)
V = V * args.scale
V = V + args.translate

x0, x1, x2 = symbols('x0 x1 x2')
m1 = Matrix([[x0, x1, x2], [1, -1, -1], [-1,1,-1], [-1,-1,1]])


vol_m1 = lambdify((x0, x1, x2), abs(m1.det()) / 6)

