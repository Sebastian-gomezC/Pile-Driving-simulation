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
flag =sys.argv[1]
parameters['allow_extrapolation'] = True
def create_mesh(mesh, cell_type, prune_z=False):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
        if prune_z:
            out_mesh.prune_z_0()
        return out_mesh
print('creando malla..')

nombre = 'zapata_refinada_prueba'


#%%
if flag == "W":
    geo_file_content= """SetFactory("OpenCASCADE");
ancho = 10;
prof =-3;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Rectangle(2) = {ancho/2-0.2, 0, 0,0.4, -0.4, 0};
Rectangle(3) = {ancho/2-0.5, -0.4, 0,1, -0.2, 0};
BooleanDifference{ Surface{1}; Delete; }{ Surface{3}; Surface{2}; Delete; }
Rectangle(4) = {0, 0, 0,ancho/2-0.5, -0.4, 0};
Rectangle(5) = {0, -0.4, 0,ancho/2-0.5, -0.2, 0};
Rectangle(6) = {0, -0.6, 0,ancho/2-0.5, prof+0.6, 0};
Rectangle(7) = {ancho/2-0.5, -0.6, 0,1, prof+0.6, 0};

Rectangle(8) = {ancho/2+0.5, 0, 0,ancho/2-0.5, -0.4, 0};
Rectangle(9) = {ancho/2+0.5, -0.4, 0,ancho/2-0.5, -0.2, 0};
Rectangle(10) = {ancho/2+0.5, -0.6, 0,ancho/2-0.5, prof+0.6, 0};
BooleanFragments{ Surface{1}; Delete; }{ Surface{4}; Surface{5}; Surface{6}; Surface{7};Surface{10}; Surface{9}; Surface{8}; Delete; }
Physical Surface("soil1",1)={4,5,6,7,8,9,10,11,12};
Physical Line("disp",1) = {2,1,13,7,8,26,27,24,22,20,15,12};
Physical Line("load",2) = {18};
Physical Line("level",3) = {19,28};
Physical Line("far",5) = {9,16,21};
//
//Transfinite Curve{8,18,12,10,17,6,16,19} = 50 ;
Transfinite Curve{28,25,23,11,9,14} = 90 Using Progression 1/1;
Transfinite Curve{19,21} = 90 Using Progression 1;
Transfinite Curve {10,12,17,22} = 48 Using Progression 1/1;
Transfinite Curve {12} = 48 Using Progression 1;
Transfinite Curve {18,16} = 20;
Transfinite Curve {15,13,26,24} = 4;
Transfinite Curve {20,4,2,7,6,27,5} = 8;
Transfinite Curve {3,1,6,8} = 6;
Transfinite Surface{4};
Transfinite Surface{5};
Transfinite Surface{6};
Transfinite Surface{7};
Transfinite Surface{8};
Transfinite Surface{9};
Transfinite Surface{10};
Transfinite Surface{11};
Transfinite Surface{12};
Mesh 2 ;
Mesh.MshFileVersion = 2.2;"""
    geo_file_content_unstruc="""SetFactory("OpenCASCADE");
ancho = 10;
prof =-3;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Rectangle(2) = {ancho/2-0.2, 0, 0,0.4, -0.4, 0};
Rectangle(3) = {ancho/2-0.5, -0.4, 0,1, -0.2, 0};
BooleanDifference{ Surface{1}; Delete; }{ Surface{3}; Surface{2}; Delete; }


Physical Line("disp",1) = {13,15,12,10,8,6};
Physical Line("load",2) = {11};
Physical Line("level",3) = {19,16};
Physical Line("far",5) = {14};
Physical Surface("soil1",1)={1};
Transfinite Curve{8,18,12,10,17,6,16,19} = 25 ;
Transfinite Curve{19} = 60 Using Progression 0.98;
Transfinite Curve{16} = 60 Using Progression 1/0.98;
Transfinite Curve{11} = 50 ;
Transfinite Curve {14} = 80 Using Bump 1.5;
Transfinite Curve{15,13} = 20 ;
Mesh 2 ;
Mesh.MshFileVersion = 2.2;
    
    """
    if os.path.exists("%s.mesh"%(nombre)):
        a=os.path.exists("%s.mesh"%(nombre))
    else:
            os.mkdir("%s.mesh"%(nombre))
    with open("%s.mesh/%s.geo"%(nombre,nombre), 'w') as filed:
        filed.write(geo_file_content)
    os.system('gmsh -2 %s.mesh/%s.geo -format msh2'%(nombre,nombre))

    # msh = meshio.read("%s.mesh/%s.msh"%(nombre,nombre))
    # triangle_mesh = create_mesh(msh, "triangle")
    # line_mesh = create_mesh(msh, "line")

    # meshio.write("mesh.xdmf", triangle_mesh)
    # meshio.write("%s.mesh/%s_mesh.xdmf"%(nombre,nombre),triangle_mesh)
    # meshio.write("%s.mesh/%s_facets.xdmf"%(nombre,nombre), line_mesh)
    os.system('dolfin-convert -i gmsh {}.mesh/{}.msh {}.mesh/{}.xml'.format(nombre, nombre,nombre,nombre))
    exit()
#%%
if flag == "R":
    mesh = Mesh()

    # mvs= MeshValueCollection("size_t", mesh, 2)
    # with XDMFFile("%s.mesh/%s_mesh.xdmf"%(nombre,nombre)) as infile:
    #     infile.read(mesh)
    #     infile.read(mvs)
    # subd = cpp.mesh.MeshFunctionSizet(mesh, mvs)

    # mvc = MeshValueCollection("size_t", mesh, 1)
    # with XDMFFile("%s.mesh/%s_facets.xdmf"%(nombre,nombre)) as infile:
    #     infile.read(mvc)
    # contorno = cpp.mesh.MeshFunctionSizet(mesh, mvc)

print('malla terminada')
#definimos la malla apartirde los archivos xml
mesh=Mesh("%s.mesh/"%(nombre)+ nombre+".xml")
tol =1E-8
cell_markers =  MeshFunction("bool", mesh,mesh.topology().dim())
cell_markers.set_all(False)

class fine(SubDomain):
        def inside(self, x, on_boundary):
            return  x[1]>=-1 and x[0]>4 and x[0]<6
fine().mark(cell_markers, True)
mesh = refine(mesh, cell_markers)
contorno = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
class load(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[1]+0.6)< tol and abs(x[0]>4.5) and abs(x[0]<5.5)
class bound(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[1]+3)< tol 
class walls(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary 
class level(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[1])<tol
walls().mark(contorno, 1)
load().mark(contorno, 2)
bound().mark(contorno, 5)
level().mark(contorno, 3)
vtkfile_u = File('%s.results/u.pvd' % (nombre))
vtkfile_fs = File('%s.results/Mohr-Coulomb_Fs.pvd' % (nombre))
vtkfile_p = File('%s.results/Pressure.pvd' % (nombre))
vtkfile_flow = File('%s.results/flow.pvd' % (nombre))

#material 

    
class K(UserExpression):
    def __init__(self, subdominio, k_0, k_1,k_2,k_3,k_4, **kwargs):
        super().__init__(**kwargs)
        self.subdominio = subdominio
        self.k_0 = k_0
        self.k_1 = k_1
        self.k_2 = k_2
        self.k_3 = k_3
        self.k_4 = k_4

    def eval_cell(self, values, x, cell):
        if self.subdominio[cell.index] == 1:
            values[0] = self.k_0
        elif self.subdominio[cell.index] == 2:
            values[0] = self.k_1
        elif self.subdominio[cell.index] == 3:
            values[0] = self.k_2
        elif self.subdominio[cell.index] == 4:
            values[0] = self.k_3
        else:
            values[0]=self.k_4

class KM(UserExpression):
    def __init__(self, subdominio, k_0, k_1,k_2,k_3,k_4, **kwargs):
        super().__init__(**kwargs)
        self.subdominio = subdominio
        self.k_0 = k_0
        self.k_1 = k_1
        self.k_2 = k_2
        self.k_3 = k_3
        self.k_4 = k_4

    def eval_cell(self, values, x, cell):
        if self.subdominio[cell.index] == 1:
            values = self.k_0
        elif self.subdominio[cell.index] == 2:
            values = self.k_1
        elif self.subdominio[cell.index] == 3:
            values = self.k_2
        elif self.subdominio[cell.index] == 4:
            values = self.k_3
        else:
            values =self.k_4

#deformacion
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    
#esfuerzo 
def sigma(u):
    return lmbda*div(u)*Identity(d) + 2*mu*epsilon(u)
ele_p  = FiniteElement("P",mesh.ufl_cell(), 1) # pressure
ele_u  = VectorElement("P",mesh.ufl_cell(), 1) # solid displacement
W = MixedElement([ele_p, ele_u])
W = FunctionSpace(mesh, W)
U = TrialFunction(W)
V = TestFunction(W)
Z= FunctionSpace(mesh, 'P', 1)
Z_e= FunctionSpace(mesh, 'DG', 0)
Z_v = VectorFunctionSpace(mesh, 'P', 1)
TS = TensorFunctionSpace(mesh, "P", 1)


p, u = split(U)
q, v = split(V)
steps =10000
n=FacetNormal(mesh)#vector normal 
t=0 # tiempo inicial
Ti=2#tiempo total
delta= Ti/steps
dt=Constant((delta))
# B_s=1E-11
# B_m=1E-10
B_f=4.4E-10
r=0.1
Poro=0.3
nu=0.3#coeficiente de poisson  
flo=Constant((0,0))
E=2729000#K(subd,2729000,2717000,334000,1252000,3286000)# #modulo elasticidad
B_m=(E/(3*(1-2*nu)))**(-1)
B_s=B_m/10
alpha=(1-B_s/B_m) #biotcoef
s_coef=(alpha-Poro)*B_s +Poro*B_f
theta =18.94#K(subd,18.94,20.61,23.27,20.53,21.84) #angulos friccion interna
C=15530#K(subd,15530,10350,18650,18400,14000) #cohesion

# E=B_m**(-1)*3*(1-2*nu)#modulo elasticidad 
print('modulo elasticidad ',E)
mu = E/(2*(1+nu))#coeficientes de Lame
rho=Constant((8000)) #densidad
lmbda = E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame




d = u.geometric_dimension()
print('geo d',d)
f = Constant((0, 0))
lam=Constant((0))
ds = Measure('ds', domain=mesh, subdomain_data=contorno)
def h(p):
    x=SpatialCoordinate(mesh)
    gam =Constant((9806.65))
    g=Constant((0,-1))
    return p/gam - x[1]

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
K1=Constant(((1E-11,0),(0,1E-11)))
K2=Constant(((6.3E-9,0),(0,6.3E-9)))
K3=Constant(((8E-8,0),(0,8E-8)))
K4=Constant(((1E-6,0),(0,1E-6)))
K5=Constant(((8E-7,0),(0,8E-7)))
visco=4E-4
K=Constant(((1E-9/visco ,0),(0,1E-9/visco )))#KM(subd,K1/visco,K2/visco,K3/visco,K4/visco,K5/visco)
H=Expression(('-gam*x[1]'),gam=gam,degree=1)

bp1=DirichletBC(W.sub(0),Constant((0)),contorno,3)


bc1 = DirichletBC(W.sub(1).sub(0), Constant((0.0)),contorno,1)
bc2 = DirichletBC(W.sub(1), Constant((0.0,0.0)),contorno,5)


x_n=Expression(('0','0','0'), gam=gam,degree=3)
X_n = Function(W)
X_n=interpolate(x_n, W)
p_n, u_n = split(X_n)

X_nn=Function(W)
X_nn=interpolate(x_n, W)
p_nn, u_nn = split(X_nn)
#
#pconst=[3./2,-2,1./2,0.0]
pconst=[1,-1,0,0]
du=pconst[0]*u
du_n=pconst[1]*u_n
du_nn=pconst[2]*u_nn
du_t= du+du_n +du_nn

divu=pconst[0]*nabla_div(u)
divu_n=pconst[1]*nabla_div(u_n)
divu_nn=pconst[2]*nabla_div(u_nn)
divu_t= divu+divu_n +divu_nn

dp=1*p
dp_n=-1*p_n
dp_nn=0*p_nn
dp_t=dp+dp_n+dp_nn


ds = Measure('ds', domain=mesh, subdomain_data=contorno)
#
#T=Expression(('0','x[1] <-I  ? 0 :  -50000*I '),I=t,r=r ,degree=2)
T=Constant((0,-50000))
F1 = inner(sigma(u), epsilon(v))*dx -alpha*p*nabla_div(v)*dx\
    -inner(T, v)*ds(subdomain_id=2, domain=mesh, subdomain_data=contorno)
F2 = dt*inner(nabla_grad(q), K*nabla_grad(p))*dx \
     + alpha*divu_t*q*dx + s_coef*(dp_t)*q*dx\
    -dt*inner(flo,nabla_grad(q))*ds(subdomain_id=5,domain=mesh, subdomain_data=contorno)-dt*inner(flo,nabla_grad(q))*ds(subdomain_id=1,domain=mesh, subdomain_data=contorno)  
L_momentum =lhs(F1)
R_momentum =rhs(F1)
L_mass=lhs(F2)
R_mass=rhs(F2)
L=L_momentum+L_mass
R=R_momentum+R_mass
bcs=[bc1,bc2,bp1]
X = Function(W)
for pot in range(steps):
    solve(L==R,X,bcs)
    print('solver end')
    X_n.assign(X)
    p_n, u_n = split(X_n)
    u_=project(u_n,Z_v)
    p_=project(p_n,Z)
    if pot % 10== 0:
        print('postproces')
        s = sigma(u_)
        
        cauchy=project(s,TS)

        o1 = sigma_1(cauchy)
        
        o2 = sigma_3(cauchy)
        
        tm=(o1-o2)/2
        fail = envFalla(o1, o2,theta,C)
        fs = tm
        fs=project(fs,Z)
        flow=-K*grad(p_)
        flow=project(flow,Z_v)
        fs.rename(" mean stress", "mean stress") ;vtkfile_fs << fs
        u_.rename("displacement", "displacement") ;vtkfile_u << u_
        flow.rename("flow", "flow") ;vtkfile_flow << flow
        p_.rename("pressure", "pressure"); vtkfile_p << p_

    print('u max:',u_.vector().get_local().max(),
              'step', pot, 'of', steps,'time:',t)
    print('p max:', p_.vector().get_local().max())
    print('p min:', p_.vector().get_local().min())
    
    t=t+delta