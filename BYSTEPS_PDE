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

@author: Sebastian
"""

import os

from fenics import*
from mshr import *
from matplotlib import pyplot
from ufl import nabla_div, nabla_grad, elem_op
import numpy as np
from ufl.tensors import as_matrix

print('creando malla..')
geo_test="""
SetFactory("OpenCASCADE");
ancho = 30 ;
prof =-20;

Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Physical Curve("disp",1) = {4};

Physical Curve("far",5) = {3,2};
Physical Surface("soil1")={1};

Characteristic Length {2,3} = 1;
Characteristic Length {1,4} = 0.03;
Mesh 2 ;

Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
"""
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

Physical Curve("far",5) = {3,7,10,13,15,16};
Physical Surface("soil1",1)={2};
Physical Surface("soil2",2)={3};
Physical Surface("soil3",3)={4};
Physical Surface("soil4",4)={5};
Physical Surface("soil5",5)={6};
Characteristic Length {4,3,6,8,10,12} = 1;
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
vtkfile_u = File('%s.results1/u.pvd' % (nombre))
vtkfile_fs = File('%s.results1/Mohr-Coulomb_Fs.pvd' % (nombre))
vtkfile_p = File('%s.results1/Pressure.pvd' % (nombre))
vtkfile_h = File('%s.results1/pressure_head.pvd' % (nombre))

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



#deformacion
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    
#esfuerzo 
def sigma(u):
    return lmbda*div(u)*Identity(d) + 2*mu*epsilon(u)
V=VectorFunctionSpace(mesh, 'CG', 1)
steps =100
n=FacetNormal(mesh)#vector normal 
t=0 # tiempo inicial
Ti=5 #tiempo total
delta= (Ti-t)/steps
dt=Constant((delta))
e=2717000 #modulo elasticidad de prueba 

E=e#K(subd,1729000,2717000,334000,1252000,3286000) #modulo elasticidad
theta =19 #K(subd,18.94,20.61,23.27,20.53,21.84)
C=18650#K(subd,15530,10350,18650,18400,14000)
nu=Constant(0.2)#coeficiente de poisson  
mu = E/2/(1+nu)#coeficientes de Lame
rho=Constant((8000)) #densidad
lmbda = E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame
u = TrialFunction(V)
v = TestFunction(V)
u_n = Function(V)
u_t = Function(V)
u_ = Function(V)
d = u.geometric_dimension()
f = Constant((0, 0))

lam=Constant((0))
Dt =Constant(dt)


ds = Measure('ds', domain=mesh, subdomain_data=contorno)
T=Constant((0,0)) #traccion exterior

def h(p):
    x=SpatialCoordinate(mesh)
    gam =Constant((9806.65))
    g=Constant((0,-1))
    return p/gam - x[1]

Q=FunctionSpace(mesh, "CG", 1)
p =TrialFunction(Q)
p_ =Function(Q)
p_n =Function(Q)
q = TestFunction(Q)
d = p.geometric_dimension()
alfa=0.85 
#theta=s/(1-s)
gam =Constant((9806.65))


def envFalla(O1,O3,Theta,c):#envolvente de falla experimental
    ang=Theta#angulo de friccion interno 
    c=c#cohesión
    return abs((O1-O3)/2 *sin(ang)+c*cos(ang))
def sigma1(sigma_x, sigma_y, tau_xy): # magnitude of first principal stress
            return((sigma_x+sigma_y)/2 + sqrt(((sigma_x-sigma_y)/2)**2 + tau_xy**2))
def sigma2(sigma_x, sigma_y, tau_xy): # magnitude of second principal stress
            return((sigma_x+sigma_y)/2 - sqrt(((sigma_x-sigma_y)/2)**2 + tau_xy**2))
def sigma_1 (T):
    c=det(T)
    b =-tr(T)
    return -b/2 + sqrt(b**2-4*c)/2
def sigma_3(T):
    c=det(T)
    b =-tr(T)
    return -b/2 - sqrt(b**2-4*c)/2
      
K=Constant(((1E-6,0),(0,1E-6)))
H=Expression(('-x[1]'),gam=gam,degree=1)
bp=DirichletBC(Q,H,contorno,5)
gamma=Constant((2))#biotcoef
r=0.15
bc1 = DirichletBC(V, Constant((0, 0)),contorno,5)
W = FunctionSpace(mesh, 'P', 1)
s_coef=5
TS = TensorFunctionSpace(mesh, "CG", 1)
for pot in range(steps):
    Dis=Expression(('x[1] <-I  ? 0 : (x[1] > -I && x[1]< -I+r ? x[1]+I: x[1]>=-I+r ? r : 0)'),I=t,r=r ,degree=1)
    bc2 = DirichletBC(V.sub(0), Dis,contorno,1)
    bcu=[bc1,bc2]
    
    
    F1 = inner(sigma(u), epsilon(v))*dx \
        - inner(f, v)*dx -\
            inner(T, v)*ds(subdomain_id=4, domain=mesh, subdomain_data=contorno)\
                +gamma*p_n*nabla_div(v)*dx- dt*nabla_div(p_n*v)*dx
    a1 = lhs(F1)
    L1 = rhs(F1)
    solve(a1==L1,u_t,bcu)
    
    
    F2=dt*inner(nabla_grad(q),K*nabla_grad(p))*dx-\
        gamma*(nabla_div(u_t)-nabla_div(u_n))*q*dx
            
    a2 = lhs(F2)
    L2 = rhs(F2)
    solve(a2==L2,p_,bp)
    
    F3 = inner(sigma(u), epsilon(v))*dx \
        - inner(f, v)*dx -\
            inner(T, v)*ds(subdomain_id=4, domain=mesh, subdomain_data=contorno)\
                +gamma*p_*nabla_div(v)*dx- dt*nabla_div(p_*v)*dx 
    a3 = lhs(F3)
    L3 = rhs(F3)
    solve(a3==L3,u_,bcu)

    
    u_n.assign(u_)
    p_n.assign(p_)
    if pot % 10 == 0:

        s = sigma(u_)
        
        cauchy=project(s,TS)

        o1 = sigma_1(cauchy)
        
        o2 = sigma_3(cauchy)
        
        tm=(o1-o2)/2
        fail = envFalla(o1, o2,theta,C)
        fs = fail/tm
        fs=project(fs,W)
        fs.rename(" Mohr-Coulomb Fs", "Mohr-Coulomb Fs") ;vtkfile_fs << fs
        u_.rename("displacement", "displacement") ;vtkfile_u << u_
        p_.rename("pressure", "pressure"); vtkfile_p << p_ 
        colum=h(p_)
        colum=project(colum, W)
        colum.rename("pressure_head", "pressure"); vtkfile_h << colum     
    print('u max:', u_.vector().get_local().max(),
              'step', pot, 'of', steps,'time:',t)
    print('p max:', p_.vector().get_local().max())
    print('p min:', p_.vector().get_local().min())
    t=t+delta
