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
from ufl import nabla_div, nabla_grad, elem_op
    

mesh= UnitIntervalMesh(100)
ele_p  = FiniteElement("P",  mesh.ufl_cell(), 1) # pressure
ele_u  = FiniteElement("P",  mesh.ufl_cell(), 2) # solid displacement
W = MixedElement([ele_p, ele_u])
W = FunctionSpace(mesh, W)
U = TrialFunction(W)
V = TestFunction(W)
p, u = split(U)
q, v = split(V)

V = FunctionSpace(mesh, "CG", 1)
U_n2 = Function(W)
p_n2, u_n2 = split(U_n2)
ul= Constant(0.001)
def u0_boundaryL(x, on_boundary):
    tol = 1e-14
    return on_boundary and abs(x[0]) < tol
def u0_boundaryR(x, on_boundary):
    tol = 1-1e-14
    return on_boundary and abs(x[0]) > tol

class Right(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1-1e-14
        return on_boundary and abs(x[0]) > tol



Gamma_0 = Right()

# Create mesh function over cell facets
exterior_facet_domains = MeshFunction("size_t", mesh,dim=1)
exterior_facet_domains.set_all(1)
Gamma_0.mark(exterior_facet_domains, 0)




vtkfile_u = File('%s.results_1d/u.pvd' % ('elasticity'))
vtkfile_p = File('%s.results_1d/Pressure.pvd' % ('elasticity'))

bc0 = DirichletBC(W.sub(1),Constant(0), u0_boundaryL)
bc1 = DirichletBC(W.sub(0),Constant(0), u0_boundaryL)
bc2 = DirichletBC(W.sub(1),Constant(-0.1), u0_boundaryR)
steps =1000
n=FacetNormal(mesh)#vector normal 
t=0 # tiempo inicial
Ti=1 #tiempo total
delta= (Ti-t)/steps
dt=Constant((delta))
K=1
fk= Constant(10**11)
E= Constant(10**11)
A= Constant(10**(-4))
L= Constant(1)
f = Constant(0)
#u0=Expression(('0','-x[0]*x[0]'),degree=2)
#U_n = interpolate(u0,W)
U_n= Function(W)
p_n, u_n = split(U_n)
U_n2= Function(W)
#U_n2 =interpolate(u0,W)
p_n2, u_n2 = split(U_n2)
T=Constant((-1))
#ds = ds(subdomain_data= exterior_facet_domains)
F = 1*dt*inner(grad(u), grad(v))*dx- dt*p.dx(0)*v*dx + (u-u_n)*v*dx - v*f*fk*dx - dt*T*v*ds
F2= 1*dt*inner(grad(q), grad(p))*dx -(u.dx(0)-u_n.dx(0))*q*dx
I1=lhs(F)
D1=rhs(F)
I2=lhs(F2)
D2=rhs(F2)
I=I1+I2
D=D1+D2
U = Function(W)
bc=[bc1,bc0]
for i in range(steps):
    if i==0:
        plot(U_n.sub(0))
        plot(U_n.sub(1))
        plt.ylim([-1.1, 1.1])
        plt.savefig('plots/step{}.png'.format(i))
        plt.clf()
    solve(I==D,U,bc)
    U_n2.assign(U_n)
    p_n2, u_n2 = split(U_n2) 
    U_n.assign(U)
    p_n, u_n = split(U_n)
    print(i)
    if i%5 ==0:
        plot(U.sub(0))
        plot(U.sub(1))
        plt.ylim([-1.1, 1.1])
        plt.savefig('plots/step{}.png'.format(i+1))
        plt.clf()
        u_=project(u_n,V)
        p_=project(p_n,V)
        u_.rename("displacement", "displacement") ;vtkfile_u << u_
        p_.rename("pressure", "pressure"); vtkfile_p << p_

