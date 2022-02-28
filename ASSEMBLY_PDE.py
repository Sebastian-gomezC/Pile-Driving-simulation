#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 21:03:08 2021

@author: jsgomez
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 14:26:23 2021

@author: SebastianG
"""

from ast import Constant
import os

from fenics import*
from mshr import *
from matplotlib import pyplot
from ufl import nabla_div, nabla_grad, elem_op
import numpy as np
from ufl.tensors import as_matrix

print('creando malla..')
geo_file_content= """
SetFactory("OpenCASCADE");
ancho = 30 ;
prof =-50;
soil1 = -4.7;
soil2= -2.7;
soil3= -11.1;
soil4= -7.7;
soil5= -23.8;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};


Rectangle(2) = {0, 0, 0, ancho, soil1, 0};

Rectangle(3) = {0, soil1, 0, ancho, soil2, 0};

Rectangle(4) = {0, soil1+soil2, 0, ancho, soil3, 0};

Rectangle(5) = {0, soil1+soil2+soil3, 0, ancho, soil4, 0};

Rectangle(6) = {0, soil1+soil2+soil3+soil4, 0, ancho, soil5, 0};

BooleanFragments{ Surface{1}; Delete; }{ Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Delete; }
Physical Curve("disp",1) = {1,5,8,11,14};
Physical Curve("level",2) = {4};
Physical Curve("far",5) = {3,7,10,13,15,16};
Physical Surface("soil1",1)={2};
Physical Surface("soil2",2)={3};
Physical Surface("soil3",3)={4};
Physical Surface("soil4",4)={5};
Physical Surface("soil5",5)={6};
Characteristic Length {4,3,6,8,10,12} = 3;
Characteristic Length {2,1,5,7,9,11} = 0.1;
Mesh 2 ;

Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
"""

nombre = 'pile_install'

if os.path.exists("%s.mesh"%(nombre)):
        a=os.path.exists("%s.mesh"%(nombre))
else:
        os.mkdir("%s.mesh"%(nombre))
with open("%s.mesh/%s.geo"%(nombre,nombre), 'w') as filed:
    filed.write(geo_file_content)
os.system('gmsh -2 {}.mesh/{}.geo -format msh2'.format(nombre,nombre))
os.system('dolfin-convert -i gmsh {}.mesh/{}.msh {}.mesh/{}.xml'.format(nombre, nombre,nombre,nombre))
print('malla terminada')
#definimos la malla apartirde los archivos xml
mesh=Mesh("%s.mesh/"%(nombre)+ nombre+".xml")
contorno = MeshFunction("size_t", mesh,"%s.mesh/"%(nombre)+nombre+"_facet_region.xml")
subd= MeshFunction("size_t",mesh,"%s.mesh/"%(nombre)+nombre+"_physical_region.xml")
vtkfile_u = File('%s.results3/u.pvd' % (nombre))
vtkfile_fs = File('%s.results3/Mohr-Coulomb_Fs.pvd' % (nombre))
vtkfile_p = File('%s.results3/Pressure.pvd' % (nombre))
vtkfile_h = File('%s.results3/pressure_head.pvd' % (nombre))

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
ele_p  = FiniteElement("DG",  mesh.ufl_cell(), 1) # pressure
ele_u  = VectorElement("CG",  mesh.ufl_cell(), 2) # solid displacement
W = MixedElement([ele_p, ele_u])
W = FunctionSpace(mesh, W)
U = TrialFunction(W)
V = TestFunction(W)
Z= FunctionSpace(mesh, 'CG', 1)
Z_v = VectorFunctionSpace(mesh, 'CG', 1)
TS = TensorFunctionSpace(mesh, "CG", 1)
X_n = Function(W)

p, u = split(U)
q, v = split(V)
p_n, u_n = split(X_n)
X_n2=Function(W)
p_n2,u_n2=split(X_n2)

steps =1000
n=FacetNormal(mesh)#vector normal 
t=0 # tiempo inicial
Ti=20 #tiempo total
delta= (Ti-t)/steps
dt=Constant((delta))
e=2717000 #modulo elasticidad de prueba 

E=K(subd,1729000,2717000,334000,1252000,3286000) #modulo elasticidad
theta =K(subd,18.94,20.61,23.27,20.53,21.84) #angulos friccion interna
C=K(subd,15530,10350,18650,18400,14000) #cohesion
nu=Constant(0.2)#coeficiente de poisson  
mu = E/2/(1+nu)#coeficientes de Lame
rho=Constant((8000)) #densidad
lmbda = E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame

d = u.geometric_dimension()
f = Constant((0, 0))
lam=Constant((0))
ds = Measure('ds', domain=mesh, subdomain_data=contorno)
T=Constant((0,0))

def h(p):
    x=SpatialCoordinate(mesh)
    gam =Constant((9806.65))
    g=Constant((0,-1))
    return p/gam - x[1]


alfa=10000 
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
      
K1=Constant(((7E-4,0),(0,7E-4)))
K2=Constant(((6.3E-4,0),(0,6.3E-4)))
K3=Constant(((8E-3,0),(0,8E-3)))
K4=Constant(((1E-6,0),(0,1E-6)))
K5=Constant(((1E-7,0),(0,1E-7)))
K=Constant(((1E-7,0),(0,1E-7)))#KM(subd,K1,K2,K3,K4,K5)
H=Expression(('0'),gam=gam,degree=1)
bp=DirichletBC(W.sub(0),H,contorno,2)
gamma=Constant((1))#biotcoef
r=0.15
bc1 = DirichletBC(W.sub(1), Constant((0, 0)),contorno,5)
flo=Constant((0))
s_coef=1.3*0.5E-9-(gamma-1.3)*0.3


F1 = inner(sigma(u), epsilon(v))*dx \
    - inner(f, v)*dx -\
    inner(T, v)*ds(subdomain_id=1, domain=mesh, subdomain_data=contorno)\
    - gamma*p*nabla_div(v)*dx 
    
F2 = dt*inner(nabla_grad(q), K*nabla_grad(p))*dx -\
    gamma*(nabla_div(u)-nabla_div(u_n))*q*dx+s_coef*(p-p_n)*q*dx - flo*q*ds(subdomain_id=1,domain=mesh, subdomain_data=contorno)


L_momentum =lhs(F1)
R_momentum =rhs(F1)
L_mass=lhs(F2)
R_mass=rhs(F2)

L= L_mass+L_momentum
R= R_mass +R_momentum

for pot in range(steps):
    if pot != 0:
        p_n2, u_n2 = X.split(deepcopy=True)
    #T=Expression(('0','x[1] <-I  ? 0 : (x[1] > -I && x[1]< -I+r ? 1000: x[1]>=-I+r ? 1000 :0)'),I=t,r=r ,degree=2)
    Dis=Expression(('x[1] <-I  ? 0 : (x[1] > -I && x[1]< -I+r ? x[1]+I: x[1]>=-I+r ? r :0)'),I=t,r=r ,degree=1)
    bc2 = DirichletBC(W.sub(1).sub(0), Dis,contorno,1)
    bcs=[bc1,bc2,bp]
    A=assemble(L)
    b=assemble(R)
    [bc.apply(A) for bc in bcs]
    [bc.apply(b) for bc in bcs]
    X = Function(W)
    solve(A, X.vector(), b,'lu')
    #solve(L==R,X,bcs)

    p_n, u_n = X.split(deepcopy=True)
    u_=as_vector((X[1],X[0]))
    u_=project(u_,Z_v)
    p_=project(X[2],Z)
    
    if pot % 100== 0:
        
        s = sigma(u_)
        
        cauchy=project(s,TS)

        o1 = sigma_1(cauchy)
        
        o2 = sigma_3(cauchy)
        
        tm=(o1-o2)/2
        fail = envFalla(o1, o2,theta,C)
        fs = fail/tm
        fs=project(fs,Z)
        fs.rename(" Mohr-Coulomb Fs", "Mohr-Coulomb Fs") ;vtkfile_fs << fs
        u_.rename("displacement", "displacement") ;vtkfile_u << u_
        p_.rename("pressure", "pressure"); vtkfile_p << p_
        colum=h(p_)
        colum=project(colum, Z)
        colum.rename("pressure_head", "pressure"); vtkfile_h << colum     
    print('u max:', u_.vector().get_local().max(),
              'step', pot, 'of', steps,'time:',t)
    print('p max:', p_.vector().get_local().max())
    print('p min:', p_.vector().get_local().min())
    t=t+delta
