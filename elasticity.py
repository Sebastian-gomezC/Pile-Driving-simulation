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
    
t=10
u0=Expression('x[0]',degree=1)
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


bc0 = DirichletBC(V,Constant(0), u0_boundary)
steps =1000
n=FacetNormal(mesh)#vector normal 
t=0 # tiempo inicial
Ti=100 #tiempo total
delta= (Ti-t)/steps
dt=Constant((delta))

fk= Constant(10**11)
E= Constant(10**11)
A= Constant(10**(-4))
L= Constant(1)
f = Constant(0)
u_n = interpolate(u0,V)
u_n2 =interpolate(u0,V)
F = dt*dt*inner(grad(u), grad(v))*10*dx + dt*1*(u-u_n)*v*dx + 100*(u-2*u_n+u_n2)*v*dx - v*f*fk*dx
I=lhs(F)
D=rhs(F)

u = Function(V)
for i in range(steps):
    solve(I==D,u,bc0)
    u_n2.assign(u_n)
    u_n.assign(u)
    print(i)
    if i%5 ==0:
        plot(u)
        plt.ylim([-1, 1])
        plt.savefig('plots/step{}.png'.format(i))
        plt.clf()

