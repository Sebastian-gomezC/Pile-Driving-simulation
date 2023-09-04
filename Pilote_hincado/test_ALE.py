#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 14:26:23 2021
@author: SebastianG
"""

from cProfile import label
import os
from pyclbr import Function
from turtle import color
from xml.etree.ElementTree import PI
import math 
from sympy import apart
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
from dolfin import *
from matplotlib import pyplot
from ufl import nabla_div, nabla_grad, elem_op
import numpy as np
from ufl.tensors import as_matrix
import  sys
import meshio
import mshr

parameters['allow_extrapolation'] = True

def create_mesh(mesh, cell_type, prune_z=False):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
        if prune_z:
            out_mesh.prune_z_0()
        return out_mesh
print('creando malla..')

nombre = 'radial_compresion'


   


#%%

mesh=RectangleMesh(Point(0.0, 0.0), Point(1,-0.5), 30, 15,"crossed")
# cell_markers =  MeshFunction("bool", mesh,mesh.topology().dim())
# cell_markers.set_all(False)
# class fine0(SubDomain):
#         def inside(self, x, on_boundary):
#             return  (x[1]>-13) and (x[0]<=3)
# fine0().mark(cell_markers, True)
# mesh = refine(mesh, cell_markers)

# cell_markers =  MeshFunction("bool", mesh,mesh.topology().dim())
# cell_markers.set_all(False)
# class fine1(SubDomain):
#         def inside(self, x, on_boundary):
#             return  (x[0]<=2 and x[1]>=-12) 
# fine1().mark(cell_markers, True)
# mesh = refine(mesh, cell_markers)
# cell_markers2 =  MeshFunction("bool", mesh,mesh.topology().dim())
# cell_markers2.set_all(False)

# class fine2(SubDomain):
#         def inside(self, x, on_boundary):
#             return  (x[0]<=1 and x[1]>=-11) 
# fine2().mark(cell_markers2, True)
# mesh = refine(mesh, cell_markers2)
# cell_markers3 =  MeshFunction("bool", mesh,mesh.topology().dim())
# cell_markers3.set_all(False)
# class fine3(SubDomain):
#         def inside(self, x, on_boundary):
#             return (x[0]<=0.5 and x[1]>=-10) 
# fine3().mark(cell_markers3, True)
# mesh = refine(mesh, cell_markers3)

# cell_markers4 =  MeshFunction("bool", mesh,mesh.topology().dim())
# cell_markers4.set_all(False)
# class fine4(SubDomain):
#         def inside(self, x, on_boundary):
#             return (x[0]<=0.2 and x[1]>=-9) 
# fine4().mark(cell_markers4, True)
# mesh = refine(mesh, cell_markers4)

contorno = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
tol=1E-6
class L(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (abs(x[0])< tol )#and (x[1]>-I)
class R(SubDomain):
        def inside(self,x,on_boundary):
             return on_boundary and  abs(x[0]-1)< tol        

class walls(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

walls().mark(contorno, 3)
L().mark(contorno, 1)
R().mark(contorno, 2)

# simetry().mark(contorno,6)

# subdominio = MeshFunction("size_t", mesh, mesh.topology().dim())
# tol=1E-2
# class S1(SubDomain):
#         def inside(self, x, on_boundary):
#             return  (x[1])>= -3
# class S2(SubDomain):
#         def inside(self, x, on_boundary):
#             return  ((x[1])<=-3) and ((x[1])>=(-8)) 
# class S3(SubDomain):
#         def inside(self, x, on_boundary):
#             return (x[1])<=(-8)
# S2().mark(subdominio, 2)
# S1().mark(subdominio, 1)
# S3().mark(subdominio, 3)
vtkfile_u = XDMFFile("%s.results/displacement.xdmf" % (nombre))
vtkfile_fs = XDMFFile("%s.results/Mohr-Coulomb_Fs.xdmf" % (nombre))
vtkfile_p = XDMFFile("%s.results/Pressure.xdmf" % (nombre))
vtkfile_flow = XDMFFile("%s.results/flow.xdmf" % (nombre))
vtkfile_bounds = XDMFFile("%s.results/bounds.xdmf" % (nombre))
#material 

vtkfile_u.parameters["flush_output"] = True
vtkfile_fs.parameters["flush_output"] = True
vtkfile_p.parameters["flush_output"] = True
vtkfile_flow.parameters["flush_output"] = True
vtkfile_bounds.parameters["flush_output"] = True
vtkfile_u.parameters["rewrite_function_mesh"] = False
vtkfile_fs.parameters["rewrite_function_mesh"] = False
vtkfile_p.parameters["rewrite_function_mesh"] = False
vtkfile_flow.parameters["rewrite_function_mesh"] = False
vtkfile_bounds.parameters["rewrite_function_mesh"] = False
    
class K(UserExpression):
    def __init__(self, subdominio, k_0, k_1,k_2, **kwargs):
        super().__init__(**kwargs)
        self.subdominio = subdominio
        self.k_0 = k_0
        self.k_1 = k_1
        self.k_2 = k_2

    def eval_cell(self, values, x, cell):
        if self.subdominio[cell.index] == 1:
            values[0] = self.k_0
        elif self.subdominio[cell.index] == 2:
            values[0] = self.k_1
        else:
            values[0]=self.k_2

class KM(UserExpression):
    def __init__(self, subdominio, k_0, k_1,k_2, **kwargs):
        super().__init__(**kwargs)
        self.subdominio = subdominio
        self.k_0 = k_0
        self.k_1 = k_1
        self.k_2 = k_2

    def value_shape(self):
        return (2,2)

    def eval_cell(self, values, x, cell):
        if self.subdominio[cell.index] == 1:
            values = self.k_0
        elif self.subdominio[cell.index] == 2:
            values = self.k_1
        else:
            values=self.k_2
def cildiv(v):
    return Dx(v[0],0)+v[0]/x[0]+Dx(v[1],1)
#deformacion
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    
#esfuerzo 
def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(d) + 2*mu*epsilon(u)
ele_p  = FiniteElement("P",mesh.ufl_cell(), 1) # pressure
ele_Q  = VectorElement("P",mesh.ufl_cell(), 1) # pressure
ele_u  = VectorElement("P",mesh.ufl_cell(), 1) # solid displacement
W = MixedElement([ele_p, ele_u,ele_Q])
W = FunctionSpace(mesh, W)
U = TrialFunction(W)
V = TestFunction(W)
Z= FunctionSpace(mesh, 'P', 1)
Z_v = VectorFunctionSpace(mesh, 'P', 1)
TS = TensorFunctionSpace(mesh, "DG", 0)

p, u,Q = split(U)
q, v,h = split(V)
steps =180
n=FacetNormal(mesh)#vector normal 
x = SpatialCoordinate(mesh)

Ti=900#tiempo total
delta= Ti/steps
dt=Constant((delta))
t=0 # tiempo inicial
# B_s=1E-11
# B_m=1E-10
B_f=4.4E-10
Poro=0.3
r=0.4
nu=0.3#coeficiente de poisson  
flo=Constant((0,0))
E=20700E3#K(subdominio,20700E3,1E6,20720E3)# #modulo elasticidad
mu = E/(2*(1+nu))#coeficientes de Lame
lmbda = E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame
B_m=(E/(3*(1-2*nu)))**(-1)
B_s=0
alpha=(1-B_s/B_m) #biotcoef
s_coef=(alpha-Poro)*B_s +Poro*B_f
theta =18.94#K(subdominio,18.94,20.61,23.27) #angulos friccion interna
C=15530#K(subdominio,15530,10350,18650) #cohesion
k=1E-8#K(subdominio,1E-8,1E-1,1E-8)

# E=B_m**(-1)*3*(1-2*nu)#modulo elasticidad 
mu = E/(2*(1+nu))#coeficientes de Lame
lmbda = E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame




d = u.geometric_dimension()
f = Constant((0, 0))
lam=Constant((0))
ds = Measure('ds', domain=mesh, subdomain_data=contorno)


#theta=s/(1-s)
gam =Constant((9806.65))


def envFalla(O1,O3,Theta,c):#envolvente de falla experimental
    ang=Theta#angulo de friccion interno 
    c=c#cohesión
    return abs((O1-O3)/2 *sin(ang)+c*cos(ang))

def sigma_1 (T):
    c=det(T)
    b =-tr(T)
    return -b/2 + sqrt(b**2-4*c)/2
def sigma_3(T):
    c=det(T)
    b =-tr(T)
    return -b/2 - sqrt(b**2-4*c)/2


vtkfile_sudint= File('%s.results/subd.pvd' % (nombre))
vtkfile_sudint<< contorno


bp1=DirichletBC(W.sub(0),Constant((0)),contorno,2)
bq1=DirichletBC(W.sub(2).sub(0),Constant((0)),contorno,1)
bq2=DirichletBC(W.sub(2).sub(1),Constant((0)),contorno,3)
bc1 = DirichletBC(W.sub(1), Constant((0.0,0.0)),contorno,2)
bc2 = DirichletBC(W.sub(1).sub(1), Constant((0.0)),contorno,3)


x_n=Expression(('0','0','0','0','0'), gam=gam,degree=5)
X_n = Function(W)
X_n=interpolate(x_n, W)
p_n, u_n,Q_n = split(X_n)

X_nn=Function(W)
X_nn=interpolate(x_n, W)
p_nn, u_nn,Q_nn = split(X_nn)
X_nnn=Function(W)
X_nnn=interpolate(x_n, W)
p_nnn, u_nnn,Q_nnn = split(X_nnn)
#
#pconst=[3./2,-2,1./2,0.0] #bdf2
#pconst = [0.48*11/6+0.52*3/2,0.48*-3+0.52*-2,0.48*3/2+0.52*1/2,0.48*-1/3] #bdf2 op
#pconst= [11/6,-3,3/2,-1/3] #bdf 3
pconst=[1,-1,0,0] #bdf1
du=pconst[0]*u
du_n=pconst[1]*u_n
du_nn=pconst[2]*u_nn
du_nnn=pconst[3]*u_nnn
du_t= du+du_n +du_nn +du_nnn

# divu=    pconst[0]*(nabla_div(u)*x[0]    + u[0])
# divu_n=  pconst[1]*(nabla_div(u_n)*x[0]  + u_n[0])
# divu_nn= pconst[2]*(nabla_div(u_nn)*x[0] + u_nn[0])
# divu_nnn=pconst[3]*(nabla_div(u_nnn)*x[0]+ u_nnn[0])
divu=    pconst[0]*nabla_div(u)
divu_n=  pconst[1]*nabla_div(u_n)
divu_nn= pconst[2]*nabla_div(u_nn)
divu_nnn=pconst[3]*nabla_div(u_nnn)
divu_t= divu+divu_n +divu_nn+divu_nnn

dp=pconst[0]*p
dp_n=pconst[1]*p_n
dp_nn=pconst[2]*p_nn
dp_nnn=pconst[3]*p_nnn
dp_t=dp+dp_n+dp_nn+dp_nnn

ds = Measure('ds', domain=mesh, subdomain_data=contorno)
#
T=Constant((0,0))
F2 = dt*inner(nabla_grad(q),nabla_grad(p))*dx \
    + alpha*divu_t*q*dx + (s_coef)*(dp_t)*q*dx \
    -dt*(inner(Constant((0,0)),n))*q*ds(subdomain_id=(3,1),domain=mesh, subdomain_data=contorno)


F1 = inner(sigma(u), epsilon(v))*dx - alpha*p*nabla_div(v)*dx\
   -inner(T, v)*ds(subdomain_id=1, domain=mesh, subdomain_data=contorno)


# T=Constant((0,0))
# F1 =  inner(sigma(u), epsilon(v))*x[0]*dx +v[0]*lmbda*nabla_div(u)*dx + (lmbda+2*mu)*v[0]*u[0]/x[0]*dx  -alpha*p*nabla_div(v)*x[0]*dx-alpha*p*v[0]*dx\
#     -inner(T, v)*ds(subdomain_id=1, domain=mesh, subdomain_data=contorno)

F3=k*inner(nabla_grad(p),h)*dx+ inner(Q,h)*dx


# F2 = dt*k*inner(nabla_grad(q),Q)*x[0]*dx\
# + alpha*divu_t*q*dx +Constant((s_coef))*(dp_t)*q*x[0]*dx\



X = Function(W)
L_momentum =lhs(F1)
R_momentum =rhs(F1)
L_mass=lhs(F2)
R_mass=rhs(F2)
L_darcy=lhs(F3)
R_darcy=rhs(F3)
L=L_momentum+L_mass+L_darcy
R=R_momentum+R_mass+R_darcy

for pot in range(steps):
    
    disp = Expression(('t/10000'),t=t,degree=1)
    bc3 = DirichletBC(W.sub(1).sub(0), disp,contorno,1)
    bcs=[bc1,bc2,bc3,bp1,bq1,bq1]
    A=assemble(L)
    b=assemble(R)
    [bc.apply(A) for bc in bcs]
    [bc.apply(b) for bc in bcs]
    print('solver')
    solve(A, X.vector(), b)
    print('solver end')
    X_nnn.assign(X_nn)
    p_nnn, u_nnn,Q_nnn = split(X_nnn)
    X_nn.assign(X_n)
    p_nn, u_nn,Q_nn = split(X_nn)
    X_n.assign(X)
    p_n, u_n,Q_n = split(X_n)
    u_=project(u_n,Z_v)
    Q_=project(Q_n,Z_v)
    p_=project(p_n,Z)
    if pot % 1== 0:
        print('postproces')
        s = sigma(u_)
        
        cauchy=project(s,TS)

        o1 = sigma_1(cauchy)
        
        o2 = sigma_3(cauchy)
        
        tm=(o1-o2)/2
        fail = envFalla(o1, o2,theta,C)
        fs = tm
        fs=project(fs,Z)
        flow=-k*Q_
        flow=project(flow,Z_v)
        fs.rename(" mean stress", "mean stress") ;vtkfile_fs.write(fs, t)
        u_.rename("displacement", "displacement") ;vtkfile_u.write(u_, t)
        flow.rename("flow", "flow") ;vtkfile_flow.write(flow, t)
        p_.rename("pressure", "pressure"); vtkfile_p.write(p_, t)
    print('u max:',u_.vector().get_local().max(),'step', pot, 'of', steps,'time:',t)
    print('p max:', p_.vector().get_local().max())
    print('p min:', p_.vector().get_local().min())
    
    t=t+delta