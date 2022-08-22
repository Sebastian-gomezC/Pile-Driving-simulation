#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 14:26:23 2021
@author: SebastianG
"""

import os
from pyclbr import Function
from turtle import color
from xml.etree.ElementTree import PI
import math 
from sympy import apart
import matplotlib.pyplot as plt
from dolfin import *
from matplotlib import pyplot
from ufl import nabla_div, nabla_grad, elem_op
import numpy as np
from ufl.tensors import as_matrix
import  sys
import meshio
parameters['allow_extrapolation'] = True
flag = sys.argv[1]
def create_mesh(mesh, cell_type, prune_z=False):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
        if prune_z:
            out_mesh.prune_z_0()
        return out_mesh
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
Transfinite Curve{1,3} = 80 ;
Transfinite Curve {2} = 40 Using Progression 1.05;
Transfinite Curve {4} = 40 Using Progression 1/1.05;
Transfinite Surface{1};
Mesh 2 ;
Mesh.MshFileVersion = 2.2;
Save StrCat(StrPrefix(General.FileName), ".msh");
"""
nombre = 'Terzagui'

# if os.path.exists("%s.mesh"%(nombre)):
#         a=os.path.exists("%s.mesh"%(nombre))
# else:
#         os.mkdir("%s.mesh"%(nombre))
# with open("%s.mesh/%s.geo"%(nombre,nombre), 'w') as filed:
#     filed.write(test_geo)
if flag == 'W':
    msh = meshio.read("%s.mesh/%s.msh"%(nombre,nombre))
    triangle_mesh = create_mesh(msh, "triangle", True)
    line_mesh = create_mesh(msh, "line", True)
    meshio.write("%s.mesh/%s_mesh.xdmf"%(nombre,nombre),triangle_mesh)
    meshio.write("%s.mesh/%s_facets.xdmf"%(nombre,nombre), line_mesh)
elif flag == 'R':
    mesh = Mesh()
    mvc = MeshValueCollection("size_t", mesh, 2)
    with XDMFFile("%s.mesh/%s_mesh.xdmf"%(nombre,nombre)) as infile:
        infile.read(mesh)
        infile.read(mvc)
    subdominio = cpp.mesh.MeshFunctionSizet(mesh, mvc)
    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("%s.mesh/%s_facets.xdmf"%(nombre,nombre)) as infile:
        infile.read(mvc)
    contorno = cpp.mesh.MeshFunctionSizet(mesh, mvc)
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
ele_p  = FiniteElement("CG",mesh.ufl_cell(), 0) # pressure
ele_u  = VectorElement("CG",mesh.ufl_cell(), 1) # solid displacement
ele_q  = VectorElement("CG",mesh.ufl_cell(), 1) # flow
W = MixedElement([ele_p, ele_u,ele_q])
W = FunctionSpace(mesh, W)
U = TrialFunction(W)
V = TestFunction(W)
Z= FunctionSpace(mesh, 'P', 1)
Z_v = VectorFunctionSpace(mesh, 'CG', 1)
TS = TensorFunctionSpace(mesh, "CG", 1)


p, u,q = split(U)
g, v,w = split(V)
steps =5000
n=FacetNormal(mesh)#vector normal 
t=0 # tiempo inicial
Ti=0.05#tiempo total
delta= Ti/steps
dt=Constant((delta))
# B_s=1E-11
# B_m=1E-10
B_f=4.4E-10

Poro=0.05
nu=0.25#coeficiente de poisson  
flo=Constant((0,0))
#E=1729000#K(subd,1729000,2717000,334000,1252000,3286000) #modulo elasticidad
B_m=1.0E-10
B_s=1.0E-11
alpha=(1-B_s/B_m) #biotcoef
s_coef=(alpha-Poro)*B_s +Poro*B_f
theta =18.94#K(subd,18.94,20.61,23.27,20.53,21.84) #angulos friccion interna
C=15530#K(subd,15530,10350,18650,18400,14000) #cohesion

E=B_m**(-1)*3*(1-2*nu)#modulo elasticidad 
print('modulo elasticidad ',E)
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
k=1.12E-14/8.9E-4
K=Constant(((k,0),(0,k)))#KM(subd,K1,K2,K3,K4,K5)
H=Expression(('-gam*x[1]'),gam=gam,degree=1)

bp1=DirichletBC(W.sub(0),Constant((0)),contorno,2)


bc1 = DirichletBC(W.sub(1).sub(0), Constant((0.0)),contorno,1)
bc2 = DirichletBC(W.sub(1), Constant((0.0,0.0)),contorno,5)


bc2 = DirichletBC(W.sub(2), Constant((0.0,0)),contorno,1)

x_n=Expression(('5E9','0','0',"0","0"), gam=gam,degree=5)
X_n = Function(W)
X_n=interpolate(x_n, W)
p_n, u_n,q_n = split(X_n)

X_nn=Function(W)
X_nn=interpolate(x_n, W)
p_nn, u_nn,q_nn= split(X_nn)
#
pconst=[3./2,-2,1./2,0.0]
#pconst=[1,-1,0,0]
du=pconst[0]*u
du_n=pconst[1]*u_n
du_nn=pconst[2]*u_nn
du_t= du+du_n +du_nn

dp= pconst[0]*p
dp_n=pconst[1]*p_n
dp_nn=pconst[2]*p_nn
dp_t=dp+dp_n+dp_nn


ds = Measure('ds', domain=mesh, subdomain_data=contorno)
T=Constant((0,-5E9))

F1 = inner(sigma(u), epsilon(v))*dx - alpha*p*nabla_div(v)*dx\
    -inner(T, v)*ds(subdomain_id=2, domain=mesh, subdomain_data=contorno)
    #+ 1*inner((sigma(u)-sigma(u_n)),epsilon(v))*dx
F2 = dt*inner(nabla_grad(g), q)*dx \
     + alpha*(nabla_div(du_t))*g*dx + s_coef*(dp_t)*g*dx
        #-dt*inner(flo,nabla_grad(q))*ds(subdomain_id=5,domain=mesh, subdomain_data=contorno)  -dt*inner(flo,nabla_grad(q))*ds(subdomain_id=1,domain=mesh, subdomain_data=contorno) 
F3 = -inner(q,w)*dx+inner(K*nabla_grad(p),w)*dx
total=F1+F2+F3
L=lhs(total)
R=rhs(total)
snaps=200
mv=1/(1/B_m+((6*(1-2*nu))/(B_m*(1+nu))))
cv=k/(alpha**2*mv+s_coef)
p0= (alpha*mv)/(alpha**2*mv+s_coef)*5E9
def p_analitico(time,snaps,cv,p0):
    z=0
    for n in range(snaps):
        
        for i in range(400):
                if i==0:
                    p=0
                else:
                    c=math.sin((2*i-1)*math.pi*(z)/(2))
                    b=math.exp(-(((2*i-1)*math.pi)/(2*0.1))**2*cv*time)
                    d=((1)/(2*i-1))
                    p=p+d*c*b
        p=4/math.pi*(p)
        if n==0:
            a=np.array([[p,z]])
        else:
            a= np.append(a,np.array([[p,z]]),axis=0)
        z=z+1/snaps
    return a

X = Function(W)
bcs=[bc1,bc2,bp1]
f = plt.figure()
f.set_figwidth(10)
f.set_figheight(10)
for pot in range(steps):

    #A=assemble(L)
    #b=assemble(R)
    #[bc.apply(A) for bc in bcs]
    #[bc.apply(b) for bc in bcs]

    solve(L==R,X,bcs)
    X_nn.assign(X_n)
    p_nn, u_nn,q_nn = split(X_nn)
    X_n.assign(X)
    p_n, u_n, q_n = split(X_n)
    u_=project(u_n,Z_v)
    q_=project(q_n,Z_v)
    p_=project(p_n,Z)
    if pot==0:
        p_0=p_(0.00,-0.1)
    tdot=(cv/0.1**2)*t
    
    if near(tdot,0.01,0.0001) or near(tdot,0.02,0.0001) or near(tdot,0.05,0.0001) or near(tdot,0.1,0.0001) or near(tdot,0.5,0.0001)or near(tdot,1,0.0001)or near(tdot,2,0.0001) :
        s = sigma(u_)
        cauchy=project(s,TS)

        o1 = sigma_1(cauchy)
        
        o2 = sigma_3(cauchy)
        
        tm=(o1-o2)/2
        fail = envFalla(o1, o2,theta,C)
        #fs = fail/tm
        z_=0
        for k in range(snaps):
            pdot=p_(0.0,z_*0.1)/p_0
            if k ==0:
                results=np.array([[pdot,-z_]])
            else:
                results =np.append(results,np.array([[pdot,-z_]]),axis=0)
            z_=z_- 1/snaps
        plt.plot(results[:, 0], results[:, 1], "-",color='silver')
        plt.xlabel("P* ")
        plt.ylabel("z*")
        p_a = p_analitico(t,snaps,cv,p0)
        print(p_a)
        plt.plot(p_a[:, 0], p_a[:, 1], "--",color='black')
        #fs.rename(" mean stress", "mean stress") ;vtkfile_fs << fs
        u_.rename("displacement", "displacement") ;vtkfile_u << u_
        q_.rename("flow", "flow") ;vtkfile_flow << q_
        p_.rename("pressure", "pressure"); vtkfile_p << p_
    
    print('u max:',u_.vector().get_local().min(),'step', pot, 'of', steps,'time*:',tdot)
    print('p max:', p_.vector().get_local().max()/p_0)
    print('p min:', p_.vector().get_local().min())
    
    t=t+delta
plt.savefig('plot_compare.png')
plt.close()
