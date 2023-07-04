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
from mpmath import mpf, besseli, besselk
flag ="R"
parameters['allow_extrapolation'] = True
nombre = 'finite_confined_aquifer'
def create_mesh(mesh, cell_type, prune_z=False):
        cells = mesh.get_cells_type(cell_type)
        cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
        out_mesh = meshio.Mesh(points=mesh.points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
        if prune_z:
            out_mesh.prune_z_0()
        return out_mesh
print('creando malla..')
R_=300
H=10
mesh=RectangleMesh(Point(0.0, 0.0), Point(R_,H), int(R_/2), int(H/2),"crossed")



contorno = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
tol =1E-6
class Q(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[0])< tol 
class confinement(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and( (abs(x[1])<tol) or (abs(x[1]-H)<tol))
class aquifer(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (abs(x[0]-R_)<tol)
confinement().mark(contorno, 2)
Q().mark(contorno, 1)
aquifer().mark(contorno, 3)
# archivos de salida 

vtkfile_contorno = File('%s.results3/subd.pvd' % (nombre))
contorno.rename("fronteras", "fronteras") ;vtkfile_contorno << contorno

vtkfile_u = File('%s.results3/u.pvd' % (nombre))
vtkfile_fs = File('%s.results3/Mohr-Coulomb_Fs.pvd' % (nombre))
vtkfile_p = File('%s.results3/Pressure.pvd' % (nombre))
vtkfile_flow = File('%s.results3/flow.pvd' % (nombre))




#deformacion
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    
#esfuerzo 
def sigma(u):
    return lmbda*div(u)*Identity(d)+2*mu*epsilon(u)
#Elementos finitos
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

steps=100
#material 
t=0 # tiempo inicial
Ti=1.1#tiempo total
delta= Ti/steps
dt=Constant((delta))
# B_s=1E-11
# B_m=1E-10
B_f=4.4E-10
Poro=0.3
nu=1/8#coeficiente de poisson  

E=30000000# #modulo elasticidad
B_m=(E/(3*(1-2*nu)))**(-1)
B_s=B_m/10
alpha=(1-B_s/B_m) #biotcoef
s_coef=(alpha-Poro)*B_s +Poro*B_f
theta =18.94#K(subd,18.94,20.61,23.27,20.53,21.84) #angulos friccion interna
C=15530#K(subd,15530,10350,18650,18400,14000) #cohesion

#E=B_m**(-1)*3*(1-2*nu)#modulo elasticidad 
print('modulo elasticidad ',E)
mu = E/(2*(1+nu))#coeficientes de Lame
print("relacion K/G= ",1/(B_m*mu))
lmbda = E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame


K=4.5E-3

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


cv_dot =(K*(1/B_m+mu/3))/(alpha**2+s_coef*(1/B_m+mu/3))

#boundary conditions 
Caudal=10
flo=Constant((Caudal/(2*math.pi*R_*H),0))
bp1=DirichletBC(W.sub(0), 0.0, "near(x[0],298,0.5) && near(x[1],6,0.5)", method="pointwise")


bc1 = DirichletBC(W.sub(1),(0.0,0.0),contorno,1)
bc2 = DirichletBC(W.sub(1).sub(1), Constant((0.0)),contorno,2)



x_n=Expression(('std::log(x[0]+1)-std::log(R)','0','0'),R=R_,degree=3)
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
scheme=types[0]
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
T=Constant((0,0))

F1 = inner(sigma(u), epsilon(v))*dx -alpha*p*nabla_div(v)*dx\
    -inner(T, v)*ds
F2 = dt*inner(nabla_grad(q), K*nabla_grad(p))*dx \
+ alpha*divu_t*q*dx + Constant((s_coef))*(dp_t)*q*dx\
-dt*inner(flo,n)*q*ds(subdomain_id=1,domain=mesh, subdomain_data=contorno)-dt*inner(flo,n)*q*ds(subdomain_id=3,domain=mesh, subdomain_data=contorno) -dt*(inner(Constant((0,0)),n))*q*ds(subdomain_id=2,domain=mesh, subdomain_data=contorno) 
L_momentum =lhs(F1) 
R_momentum =rhs(F1)
L_mass=lhs(F2)
R_mass=rhs(F2)
L=L_momentum+L_mass
R_S=R_momentum+R_mass
snaps=100

steps=100
Betta= alpha/(alpha**2 + s_coef*(1/B_m+mu/3))
q = Caudal/(2*math.pi*R_*H*K)
M=10
r = np.linspace(0.01, R_, snaps)
def pbar(r,s,f_const):
    ro=r/R_
    sigmma = (s * R_ ** 2) / cv_dot
    I_cero_ro = besseli(0, ro * (sigmma) ** (1 / 2))
    K_cero_ro = besselk(0, ro * (sigmma) ** (1 / 2))
    I_cero = besseli(0, (sigmma) ** (1 / 2))
    K_cero = besselk(0, (sigmma) ** (1 / 2))
    term_1 = (Betta * mu * f_const / q) * (1 - I_cero_ro / I_cero)
    term_2 = (1 / sigmma) *( K_cero_ro-(K_cero / I_cero) * I_cero_ro)
    v_ray = (term_1 - term_2 )
    return v_ray

def f_coef(r,s,f_const):
    ro=r/R_
    sigmma = (s * R_ ** 2)/cv_dot
    part_1 = (2*alpha/sigmma**2)*(1-1/(besseli(0,np.sqrt(sigmma))))
    part_2 = 3/(B_m*mu) + alpha*Betta*(1-(2*besseli(1,np.sqrt(sigmma)))/(np.sqrt(sigmma)*besseli(0,np.sqrt(sigmma))))
    return (part_1/part_2)*q/mu
def rad_displacement(r,s,f_const):
    ro = r/R_
    sigmma=(s * R_ ** 2) / cv_dot
    I_uno_ro = besseli(1, ro * (sigmma) ** (1 / 2))
    K_uno_ro = besselk(1, ro * (sigmma) ** (1 / 2))
    I_cero = besseli(0, (sigmma) ** (1 / 2))
    K_cero = besselk(0, (sigmma) ** (1 / 2))
    part_1_1 = ro/4- I_uno_ro/(2*sigmma**(1/2)*I_cero)
    part_1_2 = ro/4*(1/(mu*B_m)+4/3)
    part_1 = (alpha*Betta*part_1_1-part_1_2)*(mu*f_const/q)
    part_2 = (alpha/(2*sigmma**(3/2)))*(K_uno_ro - 1/(ro*sigmma ** (1/2)) + (K_cero/I_cero)*I_uno_ro)
    n_ray = part_1 + part_2
    return n_ray

def talbot(M, t,Fun,r=0,f_const=0):
    v=0
    for k in range(M - 1):
        if k == 0:
            delta =2 * M / 5
            gamma = np.exp(delta) / 2
        else:
            cot = 1 / np.tan(k * np.pi / M)
            delta = (2 * k * np.pi / 5) * (cot + 1j)
            gamma = (1 + 1j * (k * np.pi / M) *(1+ (cot) ** 2) - 1j* cot) * np.exp(delta)
        term_fun = delta / t
        termino = (gamma * Fun(r,term_fun,f_const))
        v += termino.real
    v_final=((2) / (5 * t)) * v
    return float(v_final)
    

def p_analitico(r,M,t):
    print("f_coef init")
    f_const = talbot(M, t, f_coef)
    print(f_const)
    print("p analico init")
    for l in range(len(r)):   
        if l==0:
            v_final=np.array([[cv_dot*talbot(M, t, pbar,r=r[0],f_const=f_const)/R_**2,r[l]/R_]])
        else:
            v_final= np.append(v_final,np.array([[cv_dot*talbot(M, t, pbar,r=r[l],f_const=f_const)/R_**2,r[l]/R_]]),axis=0)
    return v_final


def u_analitico(r, M, t):
    print("u analico init")
    u_final = []
    for l in range(len(r)):
        f_const = talbot(M, t, f_coef)
        u_final.append([-talbot(M, t, rad_displacement, r=r[l], f_const=f_const), r[l]/R_])
    u_final = np.array(u_final)
    return u_final


t=delta
X_w = Function(W)
bcs=[bc1,bc2,bp1]
f = plt.figure()
f.set_figwidth(10)
f.set_figheight(10)
L2=[]
uplot=[]
ig, ax = plt.subplots()
dtdot=(cv_dot/R_**2)*delta

data_001 = []
data_05 = []
data_1 = []

radio_001 = []
radio_05 = []
radio_1 = []

for pot in range(steps):


    #A=assemble(L)
    #b=assemble(R)
    #[bc.apply(A) for bc in bcs]
    #[bc.apply(b) for bc in bcs]

    solve(L==R_S,X_w,bcs)
    X_nnn.assign(X_nn)
    p_nnn, u_nnn = split(X_nnn)
    X_nn.assign(X_n)
    p_nn, u_nn = split(X_nn)
    X_n.assign(X_w)
    p_n, u_n = split(X_n)
    u_=project(u_n,Z_v)
    p_=project(p_n,Z)
    tdot= (cv_dot/R_**2)*t
    print("tiempo adimensional = ", tdot)
    # post proceso 
    

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
    z_=0
    u_analiticos=u_analitico(r,M,tdot)
    for i, row in enumerate(u_analiticos):  
        if near(row[1], 0.01, dt / 2):  
            print("YUJU")
            radio_001.append([row[0], tdot])
        elif near(row[1], 0.2, dt / 2):
            print("uy")
            radio_05.append([row[0], tdot])
        elif near(row[1], 1, dt / 2):
            print("caray")
            radio_1.append([row[0], tdot])
         
    if near(t,0.01/(cv_dot/R_**2),dt/2):
        print(tdot)
        data_001.append(u_analiticos)
    elif near(t,0.2/(cv_dot/R_**2),dt/2):
        print(tdot)
        data_05.append(u_analiticos)
    elif near(t,0.99/(cv_dot/R_**2),dt/2):
        print(tdot)
        data_1.append(u_analiticos)

    if near(t,0.01/(cv_dot/R_**2),dt/2) or near(t,0.1/(cv_dot/R_**2),dt/2) or near(t,1/(cv_dot/R_**2),dt/2) :
        for y in range(snaps):
            pdot=p_(z_*R_,5)/q
            if y ==0:
                results=np.array([[pdot,z_]])
            else:
                results =np.append(results,np.array([[pdot,z_]]),axis=0)
            z_+=1/snaps
        p_a = p_analitico(r,M,t)
        L2.append([np.sum((results[:,0]-p_a[:,0])**2),tdot])
        
        line1, =ax.plot(results[:, 1],results[:, 0], "-",color='red',label='FEM',)
        plt.xlabel("$r/R$ ",fontsize=20)
        plt.ylabel("$p/q$",fontsize=20)
        line2,=ax.plot(p_a[:, 1], p_a[:, 0], "--",color='black',label='Analítica')
        
        lines = plt.gca().get_lines()
        l1=lines[-1]
        print("labelproblem")
        labelLine(l1,0.1,label=r'$t*=${}'.format(round(tdot,2)),ha='left',va='bottom',align = True)
        print("labelproblem")
        if near(t,1/(cv_dot/R_**2),dt/2) :
            ax.legend(handler_map={line1: HandlerLine2D(numpoints=4)},loc= 'upper center')
            
        #uplot.append([u_analitico(r,M,t),(1/B_m+mu/3)*u_(R_,0)[0]/(R_*q),tdot])
       # print('error en la norma l2',np.sum((results[:,0]-p_a[:,0])**2 ))
        
      
        
    if near(t,1/(cv_dot/R_**2),dt/2):
        break

    print('step', pot, 'of', steps,'time*:',tdot)

    t += delta
    
plt.savefig('Acuifero_confinado/RESULTADOS_ACUIFERO/presion_aquifero_dt%s_grosse_%s.png'%(round(dtdot,5),scheme),dpi=300)
plt.close()
print(p_)

 
u_final=u_analitico(r,M,tdot)
print(u_final)

data_001 = np.concatenate(data_001)
data_05 = np.concatenate(data_05)
data_1 = np.concatenate(data_1)

# Trazar la gráfica con las tres líneas
plt.plot(data_001[:, 1], data_001[:, 0], label='t*=0.01')
plt.plot(data_05[:, 1], data_05[:, 0], label='t*=0.2')
plt.plot(data_1[:, 1], data_1[:, 0], label='t*=1')
plt.xlabel("$r/R$")
plt.ylabel("$u_z$")
plt.legend()
plt.grid(True, color='k', which="both", alpha=0.3, linestyle='-', linewidth=0.5)
plt.savefig('Acuifero_confinado/RESULTADOS_ACUIFERO/desp_aquifero_dt%s_grosse_%s.png'%(round(dtdot,5),scheme),dpi=300)
plt.close()
print(radio_001)

x_values001 = [item[0] for item in radio_001]
y_values001 = [item[1] for item in radio_001]
x_values05 = [item[0] for item in radio_05]
y_values05 = [item[1] for item in radio_05]
x_values1 = [item[0] for item in radio_1]
y_values1 = [item[1] for item in radio_1]
# Trazar la gráfica con las tres líneas
plt.plot(x_values001, y_values001, label='t*=0.01')
plt.plot(x_values05, y_values05, label='t*=0.2')
plt.plot(x_values1, y_values1, label='t*=1')
plt.xlabel("$t*$")
plt.ylabel("$u_z$")
plt.legend()
plt.grid(True, color='k', which="both", alpha=0.3, linestyle='-', linewidth=0.5)
plt.savefig('Acuifero_confinado/RESULTADOS_ACUIFERO/desp(pvst)_aquifero_dt%s_grosse_%s.png'%(round(dtdot,5),scheme),dpi=300)
plt.close()