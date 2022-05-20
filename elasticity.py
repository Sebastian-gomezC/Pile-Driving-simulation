#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 0:30:26 2020
@author: sebastian
"""

from __future__ import print_function
from ast import Constant, Expression
from pyclbr import Function
from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
    
t=1

ul=0.0
f=10**11
E=10**11
A=10**(-4)
L= 1

mesh= UnitIntervalMesh(t)
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)
u_n=Function(V)
u_n2 = Function(V)
ul= Constant(0.001)
def u0_boundary(x, on_boundary):
    tol = 1e-14
    return on_boundary and abs(x[0]) < tol

def u0_boundaryR(x, on_boundary):
    tol = 1-1e-14
    return on_boundary and abs(x[0]) > tol

class Right(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1-1e-3
        return on_boundary and abs(x[0]) > tol



Gamma_0 = Right()

# Create mesh function over cell facets
exterior_facet_domains = MeshFunction("size_t", mesh,dim=1)
exterior_facet_domains.set_all(1)
Gamma_0.mark(exterior_facet_domains, 1)

bc0 = DirichletBC(V,Constant(0), u0_boundary)
steps =100
n=FacetNormal(mesh)#vector normal 
t=0 # tiempo inicial
Ti=1 #tiempo total
delta= (Ti-t)/steps
dt=Constant((delta))

fk= Constant(1)
E= Constant(1)
A= Constant(1)
L= Constant(1)
f = Constant(10)
u0=Expression('10*x[0]*x[0]',degree=1)
#u_n = interpolate(u0,V)
u_n = Function(V)
#u_n2 =interpolate(u0,V)
u_n2 = Function(V)
#ds = ds(subdomain_data= exterior_facet_domains
T=Constant((-1000))
Fs=inner(grad(u), grad(v))*60*dx + - dt*T*v*ds
Is=lhs(Fs)
Ds=rhs(Fs)
us = Function(V)
solve(Is==Ds,us,bc0)
plot(us)
plt.ylim([-0.18, 0.18])
plt.savefig('plots1/step{}.png'.format(-1))
plt.clf()
hdfw = HDF5File(mesh.mpi_comm(), "a.h5", "w")
hdfw.write(mesh, "mesh2")
hdfw.write(us, "despl")
hdfw.close()
F = dt*dt*inner(grad(u), grad(v))*60*dx + dt*10*(u-u_n)*v*dx  
I=lhs(F)
D=rhs(F)


mesh2 = Mesh()
hdf5 = HDF5File(mesh2.mpi_comm(), "a.h5", "r")
hdf5.read(mesh2, "mesh2")
hdf5.read(u_n, "despl")
hdf5.read(u_n2, "despl")
hdf5.close()
u = Function(V)
for i in range(steps):
    solve(I==D,u,bc0)
    # if i ==10:
    #     T.assign(Constant(0))
    #     I=lhs(F)
    #     D=rhs(F)
    # u_n2.assign(u_n)
    # u_n.assign(u)
    print(i)
    # if i%5 ==0:
    plot(u)
    plt.ylim([-0.18, 0.18])
    plt.savefig('plots1/step{}.png'.format(i))
    plt.clf()
