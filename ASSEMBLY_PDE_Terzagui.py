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

test_geo="""
SetFactory("OpenCASCADE");
ancho = 0.05;
prof =-0.1;
Rectangle(1) = {0, 0, 0, ancho, prof, 0};
Physical Line("disp",1) = {4,2};
Physical Line("level",2) = {1};
Physical Line("far",5) = {3};
Physical Surface("soil1",1)={1};
Transfinite Curve{1,2,3,4} = 31 ;
Transfinite Surface{1};
Mesh 2 ;
Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
"""
nombre = 'Terzagui'

if os.path.exists("%s.mesh"%(nombre)):
        a=os.path.exists("%s.mesh"%(nombre))
else:
        os.mkdir("%s.mesh"%(nombre))
with open("%s.mesh/%s.geo"%(nombre,nombre), 'w') as filed:
    filed.write(test_geo)
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
vtkfile_flow = File('%s.results3/flow.pvd' % (nombre))

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
ele_p  = FiniteElement("P",  mesh.ufl_cell(), 1) # pressure
ele_u  = VectorElement("P",  mesh.ufl_cell(), 2) # solid displacement
W = MixedElement([ele_p, ele_u])
W = FunctionSpace(mesh, W)
U = TrialFunction(W)
V = TestFunction(W)
Z= FunctionSpace(mesh, 'P', 1)
Z_v = VectorFunctionSpace(mesh, 'CG', 1)
TS = TensorFunctionSpace(mesh, "CG", 1)


p, u = split(U)
q, v = split(V)



steps =100
n=FacetNormal(mesh)#vector normal 
t=0 # tiempo inicial
Ti=1#tiempo total
delta= Ti/steps
dt=Constant((delta))
# B_s=1E-11
# B_m=1E-10
B_f=2.2E-9

r=0.45
Poro=0.5
nu=0.3#coeficiente de poisson  
flo=Constant((0,0))
E=1729000#K(subd,1729000,2717000,334000,1252000,3286000) #modulo elasticidad
B_m=(E/(3*(1-2*nu)))**(-1)
B_s=B_m/10
gamma=Constant((1-B_s/B_m)) #biotcoef
s_coef=(gamma-Poro)*B_s +Poro*B_f
theta =18.94#K(subd,18.94,20.61,23.27,20.53,21.84) #angulos friccion interna
C=15530#K(subd,15530,10350,18650,18400,14000) #cohesion

#E=B_m**(-1)*3*(1-2*nu)#modulo elasticidad de prueba 

print('modulo elasticidad ',B_m)
mu = E/2/(1+nu)#coeficientes de Lame
rho=Constant((8000)) #densidad
lmbda = E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame




d = u.geometric_dimension()
f = Constant((0, 0))
lam=Constant((0))
ds = Measure('ds', domain=mesh, subdomain_data=contorno)


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
K5=Constant(((8E-7,0),(0,8E-7)))
K=Constant(((1.12E-8,0),(0,1.12E-8)))#KM(subd,K1,K2,K3,K4,K5)
H=Expression(('-gam*x[1]'),gam=gam,degree=1)

bp1=DirichletBC(W.sub(0),Constant((0)),contorno,2)


bc1 = DirichletBC(W.sub(1).sub(0), Constant((0.0)),contorno,1)
bc2 = DirichletBC(W.sub(1).sub(1), Constant((0.0)),contorno,5)



#X_n=Expression(('0','0','0'), gam=gam,degree=3)
#X_n=interpolate(X_n, W)
X_n = Function(W)
p_n, u_n = split(X_n)
ds = Measure('ds', domain=mesh, subdomain_data=contorno)
T=Constant((0,0))

F1 = inner(sigma(u), epsilon(v))*dx - gamma*p*nabla_div(v)*dx\
    - inner(f, v)*dx \
    - inner(T, v)*ds(subdomain_id=2, domain=mesh, subdomain_data=contorno)\
     #+ 100*inner((sigma(u)-sigma(u_n)),epsilon(v))*dx
F2 = dt*inner(nabla_grad(q), K*nabla_grad(p))*dx +\
    gamma*(nabla_div(u)-nabla_div(u_n))*q*dx +s_coef*(p-p_n)*q*dx\
        -dt*inner(flo,n)*q*ds(subdomain_id=5,domain=mesh, subdomain_data=contorno)  -dt*inner(flo,n)*q*ds(subdomain_id=1,domain=mesh, subdomain_data=contorno) 


L_momentum =lhs(F1)
R_momentum =rhs(F1)
L_mass=lhs(F2)
R_mass=rhs(F2)

L= L_mass+L_momentum
R= R_mass +R_momentum
X = Function(W)
#bcs=[bc1,bc2,bp1]
for pot in range(steps):
    if t> (0.5):
          bc3 = DirichletBC(W.sub(1).sub(1), Constant((-0.5*0.02)),contorno,2)
          bcs=[bc1,bc2,bc3,bp1]
    else:
        bc3 = DirichletBC(W.sub(1).sub(1), Constant((-t*0.02)),contorno,2)
        bcs=[bc1,bc2,bc3,bp1]
    
    
    A=assemble(L)
    b=assemble(R)
    [bc.apply(A) for bc in bcs]
    [bc.apply(b) for bc in bcs]

    solve(A, X.vector(), b,'lu')
    X_n.assign(X)
    p_n, u_n = split(X_n)
    u_=project(u_n,Z_v)
    p_=project(p_n,Z)
    if pot % 1== 0:
    
        s = sigma(u_)
        
        cauchy=project(s,TS)

        o1 = sigma_1(cauchy)
        
        o2 = sigma_3(cauchy)
        
        tm=(o1-o2)/2
        fail = envFalla(o1, o2,theta,C)
        #fs = fail/tm
        fs=project(tm,Z)
        flow=-K*grad(p_)
        flow=project(flow,Z_v)
        fs.rename(" mean stress", "mean stress") ;vtkfile_fs << fs
        u_.rename("displacement", "displacement") ;vtkfile_u << u_
        flow.rename("flow", "flow") ;vtkfile_flow << flow
        p_.rename("pressure", "pressure"); vtkfile_p << p_

    print('u max:',u_.vector().get_local().min(),
              'step', pot, 'of', steps,'time:',t)
    print('p max:', p_.vector().get_local().max())
    print('p min:', p_.vector().get_local().min())
    
    t=t+delta
