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
from label_lines import *
flag ="R"
parameters['allow_extrapolation'] = True
nombre = 'Terzaghi_multi'
def create_mesh(mesh, cell_type, prune_z=False):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
        if prune_z:
            out_mesh.prune_z_0()
        return out_mesh
print('creando malla..')
L=-3
H=1.5
mesh=RectangleMesh(Point(-H/2, 0.0), Point(H/2,L), int(H*10), int(-L*20),"crossed")
h1=-1
#geo = mshr.Rectangle(Point(-0.25, 0.0), Point(0.25,-1))
#mesh = mshr.generate_mesh(geo,25)
# cell_markers =  MeshFunction("bool", mesh,mesh.topology().dim())
# cell_markers.set_all(False)
# class fine(SubDomain):
#         def inside(self, x, on_boundary):
#             return  x[1]>=-0.1
# fine().mark(cell_markers, True)
# mesh = refine(mesh, cell_markers)

contorno = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
tol =1E-3
class load(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[1])< tol 
class bound(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[1]-L)< tol 
class walls(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and ((abs(x[0]+H/2)< tol) or (abs(x[0]-H/2)< tol) )
load().mark(contorno, 2)
bound().mark(contorno, 5)
walls().mark(contorno, 1)
# archivos de salida 
subd = MeshFunction("size_t", mesh, mesh.topology().dim())
class H1(SubDomain):
    def inside(self,x,on_boundary):
        return x[1]>=h1 
class H2(SubDomain):
    def inside(self,x,on_boundary):
        return (x[1]-h1)<tol

H2().mark(subd,2)
H1().mark(subd,1)
vtkfile_subd = File('%s.results3/subd.pvd' % (nombre))
subd.rename("subdomain", "subdomain") ;vtkfile_subd << subd

vtkfile_u = File('%s.results3/u.pvd' % (nombre))
vtkfile_fs = File('%s.results3/Mohr-Coulomb_Fs.pvd' % (nombre))
vtkfile_p = File('%s.results3/Pressure.pvd' % (nombre))
vtkfile_flow = File('%s.results3/flow.pvd' % (nombre))

#material 

    
class K(UserExpression):
    def __init__(self, subdominio, k_0, k_1, **kwargs):
        super().__init__(**kwargs)
        self.subdominio = subdominio
        self.k_0 = k_0
        self.k_1 = k_1

    def eval_cell(self, values, x, cell):
        if self.subdominio[cell.index] == 1:
            values[0] = self.k_0
        elif self.subdominio[cell.index] == 2:
            values[0] = self.k_1
class KM(UserExpression):
    def __init__(self, subdominio, k_0, k_1, **kwargs):
        super().__init__(**kwargs)
        self.subdominio = subdominio
        self.k_0 = k_0
        self.k_1 = k_1

    def eval_cell(self, values, x, cell):
        if self.subdominio[cell.index] == 1:
            values = self.k_0
        elif self.subdominio[cell.index] == 2:
            values = self.k_1

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
Z_d= FunctionSpace(mesh, 'DG', 0)
Z_v = VectorFunctionSpace(mesh, 'P', 1)
TS = TensorFunctionSpace(mesh, "P", 1)


p, u = split(U)
q, v = split(V)
steps =1000

t=0 # tiempo inicial
Ti=10#tiempo total
delta= Ti/steps
dt=Constant((delta))
# B_s=1E-11
B_m=0.1E-5 #compresibilidad modulo bulk moduls
B_f=4.4E-5
Poro=0.3
nu=0.3#coeficiente de poisson  
flo=Constant((0,0))
#E=50000000# #modulo elasticidad
#B_m=(E/(3*(1-2*nu)))**(-1)
B_s=B_m/10
alpha=(1-B_s/B_m) #biotcoef
s_coef=(alpha-Poro)*B_s +Poro*B_f
theta =18.94#K(subd,18.94,20.61,23.27,20.53,21.84) #angulos friccion interna
C=15530#K(subd,15530,10350,18650,18400,14000) #cohesion

E=B_m**(-1)*3*(1-2*nu)#modulo elasticidad 
print('modulo elasticidad ',E)
mu = E/(2*(1+nu))#coeficientes de Lame
lmbda = E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame




d = u.geometric_dimension()
f = Constant((0, 0))



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

n=FacetNormal(mesh)#vector normal      
k=1.15E-5
k1=Constant(((1.15E-5,0),(0,1.15E-5)))
k2=Constant(((1.15E-6,0),(0,1.15E-6)))
K=KM(subd,k1,k2)

bp1=DirichletBC(W.sub(0),Constant((0)),contorno,2)


bc1 = DirichletBC(W.sub(1).sub(0),(0.0),contorno,1)
bc2 = DirichletBC(W.sub(1), Constant((0.0,0.0)),contorno,5)



x_n=Expression(('0','0','0'),degree=3)
X_n = Function(W)
X_n=interpolate(x_n, W)
p_n, u_n = split(X_n)

X_nn=Function(W)
X_nn=interpolate(x_n, W)
p_nn, u_nn = split(X_nn)

X_nnn=Function(W)
X_nnn=interpolate(x_n, W)
p_nnn, u_nnn = split(X_nnn)
#
#pconst=[3./2,-2,1./2,0.0] #bdf2
#pconst = [0.48*11/6+0.52*3/2,0.48*-3+0.52*-2,0.48*3/2+0.52*1/2,0.48*-1/3] #bdf2 op
#pconst= [11/6,-3,3/2,-1/3] #bdf 3
pconst=[1,-1,0,0] #bdf1
types=['BDF1','BDF2','BDF2op','BDF3']
scheme=types[3]
du=pconst[0]*u
du_n=pconst[1]*u_n
du_nn=pconst[2]*u_nn
du_nnn=pconst[3]*u_nnn
du_t= du+du_n +du_nn +du_nnn

divu=pconst[0]*nabla_div(u)
divu_n=pconst[1]*nabla_div(u_n)
divu_nn=pconst[2]*nabla_div(u_nn)
divu_nnn=pconst[3]*nabla_div(u_nnn)
divu_t= divu+divu_n +divu_nn+divu_nnn

dp=pconst[0]*p
dp_n=pconst[1]*p_n
dp_nn=pconst[2]*p_nn
dp_nnn=pconst[3]*p_nnn
dp_t=dp+dp_n+dp_nn+dp_nnn


ds = Measure('ds', domain=mesh, subdomain_data=contorno)
T=Constant((0,-5E7))

F1 = inner(sigma(u), epsilon(v))*dx -alpha*p*nabla_div(v)*dx\
    -inner(T, v)*ds(subdomain_id=2, domain=mesh, subdomain_data=contorno)
F2 = dt*inner(nabla_grad(q), K*nabla_grad(p))*dx \
     + alpha*divu_t*q*dx + s_coef*(dp_t)*q*dx\
    -dt*inner(flo,nabla_grad(q))*ds(subdomain_id=5,domain=mesh, subdomain_data=contorno)  -dt*inner(flo,nabla_grad(q))*ds(subdomain_id=1,domain=mesh, subdomain_data=contorno) 
L_momentum =lhs(F1)
R_momentum =rhs(F1)
L_mass=lhs(F2)
R_mass=rhs(F2)
L=L_momentum+L_mass
R=R_momentum+R_mass
snaps=100
mv=1/((1/B_m)+((4/3)*mu))
cv=k/(alpha**2*mv+s_coef)
p0= ((alpha*mv)/(alpha**2*mv+s_coef))*5E7
def p_analitico(time,snaps,cv,p0):
    z=0
    for n in range(snaps):
        for i in range(100):
                if i==0:
                    p=0
                else:
                    c=math.sin((2*i-1)*math.pi*(z)/(2))
                    b=math.exp(-(((2*i-1)*math.pi)/(2*1))**2*cv*time)
                    d=((1)/(2*i-1))
                    p=p+d*c*b
        p=4/math.pi*(p)
        if n==0:
            a=np.array([[p,z]])
        else:
            a= np.append(a,np.array([[p,z]]),axis=0)
        z=z+1/snaps
    return a
def u_analitico(time,snaps,cv,p0):
    c=-mv*5e7*1
    for i in range(100):
            if i==0:
                p=0
            else:
                b=math.exp(-(((2*i-1)*math.pi)/(2*1))**2*cv*time)
                d=((1)/(2*i-1)**2)
                p+=b*d
    a=8*alpha*mv*1/(math.pi**2)*p0*(p)+c
    return a

X_w = Function(W)
bcs=[bc1,bc2,bp1]
f = plt.figure()
f.set_figwidth(10)
f.set_figheight(10)
L2=[]
uplot=[]
u_f=mv*5E7
print('u final',u_f)
ig, ax = plt.subplots()
dtdot=(cv/1**2)*delta
for pot in range(steps):
    
    #A=assemble(L)
    #b=assemble(R)
    #[bc.apply(A) for bc in bcs]
    #[bc.apply(b) for bc in bcs]

    solve(L==R,X_w,bcs)
    X_nnn.assign(X_nn)
    p_nnn, u_nnn = split(X_nnn)
    X_nn.assign(X_n)
    p_nn, u_nn = split(X_nn)
    X_n.assign(X_w)
    p_n, u_n = split(X_n)
    u_=project(u_n,Z_v)
    p_=project(p_n,Z)

    # post proceso 
    if pot==0:
        u_0dot = (mv-(alpha**2*mv**2)/(alpha**2*mv+s_coef))*5e7
        p_0=5E7#p_(0.00,-0.1)
    tdot=(cv/1**2)*t
    if pot%(steps/50) ==0:
        s = sigma(u_)
        cauchy=project(s,TS)

        o1 = sigma_1(cauchy)
            
        o2 = sigma_3(cauchy)
            
        tm=-o2#(o1-o2)/2
        fail = envFalla(o1, o2,theta,C)
        fs = tm
        fs=project(fs,Z)
        flow=-K*grad(p_)
        flow=project(flow,Z_v)
        fs.rename("mean stress", "mean stress") ;vtkfile_fs << fs
        u_.rename("displacement", "displacement") ;vtkfile_u << u_
        flow.rename("flow", "flow") ;vtkfile_flow << flow
        p_.rename("pressure", "pressure"); vtkfile_p << p_
    uplot.append([(u_analitico(t, snaps, cv, p0)+u_0dot)/(u_f-u_0dot),(u_(0,0)[1]+u_0dot)/(u_f-u_0dot),tdot])
    print(tdot)
    z_=0
    for k in range(snaps):
        pdot=p_(0.0,z_*3)/p0
        if k ==0:
            results=np.array([[pdot,-z_]])
        else:
            results =np.append(results,np.array([[pdot,-z_]]),axis=0)
        z_=z_- 1/snaps
    p_a = p_analitico(t,snaps,cv,p0)
    L2.append([np.sum((results[:,0]-p_a[:,0])**2),tdot])
    if near(t,0.01/(cv/1**2),dt/2) or near(t,0.1/(cv/1**2),dt/2) or near(t,0.1/(cv/1**2),dt/2) or near(t,0.5/(cv/1**2),dt/2)or near(t,1/(cv/1**2),dt/2):
        
        
        line1, =ax.plot(results[:, 0], results[:, 1], "-",color='red',label='FEM',)
        plt.xlabel("$P*$ ",fontsize=20)
        plt.ylabel("$z*$",fontsize=20)
        line2,=ax.plot(p_a[:, 0], p_a[:, 1], "--",color='black',label='Analítica')
        
        lines = plt.gca().get_lines()
        l1=lines[-1]
        labelLine(l1,0.4,label=r'$t*=${}'.format(round(tdot,2)),ha='left',va='bottom',align = True)
        if near(t,0.01/(cv/1**2),dt/2) :
            ax.legend(handler_map={line1: HandlerLine2D(numpoints=4)},loc= 'upper center')
            
        plt.ylim((0,1.05))
        plt.xlim((0,1.05))
        print('error en la norma l2',np.sum((results[:,0]-p_a[:,0])**2 ))
        print('u max:',u_.vector().get_local().min(),'step', pot, 'of', steps,'time*:',tdot)
        print('p max:', p_.vector().get_local().max())
        print('p min:', p_.vector().get_local().min())
        
        
    if near(t,2/(cv/1**2),dt/2):
        break
    
    t += delta
plt.savefig('resultados_%s/consolidacion_dt%s_grosse_%s.png'%(nombre,round(dtdot,5),scheme),dpi=300)
plt.close()
ig, ax = plt.subplots()
uplot=np.array(uplot)
line1, =ax.plot(uplot[:,2],-uplot[:,0],"--",color='black',label='Analítica')
line2,=ax.plot(uplot[:,2],-uplot[:,1], "-",color='red',label='FEM')
ax.legend(handler_map={line1: HandlerLine2D(numpoints=4)},loc= 'upper left')
ax.set_xscale('log')
plt.xlabel("$t*$ ",fontsize=20)
plt.ylabel("$u*_{z}$",fontsize=20)
plt.grid(True,color='k',which="both",alpha=0.3, linestyle='-', linewidth=0.5)
plt.savefig('resultados_%s/disp_dt%s_grosse_%s.png'%(nombre,round(dtdot,5),scheme),dpi=300,)
plt.close()

L2=np.array(L2)
plt.semilogy(L2[:,1],L2[:,0])
plt.savefig('resultados_%s/L2norm_dt%s_grosse_%s.png'%(nombre,round(dtdot,5),scheme),dpi=300)
np.savetxt('resultados_%s/L2norm_dt%s_grosse_%s.out'%(nombre,round(dtdot,5),scheme), (L2)) 
plt.close()

