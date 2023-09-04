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
from ufl import nabla_div, nabla_grad,bessel_I,bessel_K
import numpy as np
import  sys
import meshio
import mshr
from label_lines import *
import mpmath 
def besseli(n,f):
    return bessel_I(n, f).evaluate(1,1,1,1)
def besselk(n,f):
    if n==0:
        bes = bessel_K(n, f).evaluate(1,1,1,1)
    else:
        bes = mpmath.besselk(n,f)
    return bes
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
mesh=RectangleMesh(Point(0.0, 0.0), Point(R_,H), int((int(H/2)+1)*R_/H), int(H/2)+1,"crossed")


contorno = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
tol =1E-6
class Q(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[0]-0.0)< tol 
class down(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and abs(x[1])<tol
class up(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and abs(x[1]-H)<tol
class aquifer(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (abs(x[0]-R_)<tol)
down().mark(contorno, 2)
aquifer().mark(contorno, 3)
up().mark(contorno,4)
Q().mark(contorno, 1)
# archivos de salida 

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
vtkfile_bounds.write(contorno)

#deformacion
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    
#esfuerzo 
def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(d)+2*mu*epsilon(u)
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

# B_s=1E-11
# B_m=1E-10
B_f= 1E-9
Poro=0.4

K_m=500E3
nu=0.225#0.386#coeficiente de poisson  
#E=3*K_m*(1-2*nu)#(9*K_m*G_m)/(3*K_m+G_m)# #modulo elasticidad
B_m=(K_m)**(-1)
B_s=0
alpha=(1-B_s/B_m) #biotcoef
s_coef=(alpha-Poro)*B_s +Poro*B_f

theta =18.94#K(subd,18.94,20.61,23.27,20.53,21.84) #angulos friccion interna
C=15530#K(subd,15530,10350,18650,18400,14000) #cohesion

#E=B_m**(-1)*3*(1-2*nu)#modulo elasticidad 

mu =3*K_m*(1-2*nu)/(2*(1+nu))#coeficientes de Lame
print("relacion K/G= ",1/(B_m*mu))
lmbda =(3*K_m*nu)/(1+nu) #E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame


K=1.15741E-14

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
x = SpatialCoordinate(mesh)
cv_dot =(K*(K_m+mu/3))/(alpha**2+s_coef*(K_m+mu/3))
cv =    (K*(K_m+4*mu/3))/(alpha**2+s_coef*(K_m+4*mu/3))

print('cv*=   ',cv_dot)
print('cv=    ',cv)
#material 
t=0 # tiempo inicial



dtdot=0.001
delta= dtdot/(cv_dot/R_**2)
steps=int( 2/dtdot)
dt=Constant((delta))
#boundary conditions 
print('consolidation time =',1*R_**2/cv_dot)
Caudal=1E-5
flo=Constant((Caudal/(2*math.pi*H)))
#bp1=DirichletBC(W.sub(0), 0.0, "near(x[0],298,0.5) && near(x[1],6,0.5)", method="pointwise")
bp1=DirichletBC(W.sub(0), (0.0), contorno,3)
bc1 = DirichletBC(W.sub(1).sub(0),Constant(0.0),contorno,1)
bc2 = DirichletBC(W.sub(1).sub(1), Constant((0.0)),"near(x[1],0,0.01)", method="pointwise")



x_n=Expression(('0','0','0'),R=R_,degree=3)
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
types=['BDF1_final','BDF2','BDF2op','BDF3_comp_cil']
scheme=types[0]
du=pconst[0]*u
du_n=pconst[1]*u_n
du_nn=pconst[2]*u_nn
du_nnn=pconst[3]*u_nnn
du_t= du+du_n +du_nn +du_nnn

divu=    pconst[0]*(nabla_div(u)*x[0]+ u[0])
divu_n=  pconst[1]*(nabla_div(u_n)*x[0]+ u_n[0])
divu_nn= pconst[2]*(nabla_div(u_nn)*x[0]+ u_nn[0])
divu_nnn=pconst[3]*(nabla_div(u_nnn)*x[0]+ u_nnn[0])
divu_t= divu+divu_n +divu_nn+divu_nnn

dp=pconst[0]*p
dp_n=pconst[1]*p_n
dp_nn=pconst[2]*p_nn
dp_nnn=pconst[3]*p_nnn
dp_t=dp+dp_n+dp_nn+dp_nnn


ds = Measure('ds', domain=mesh, subdomain_data=contorno)
T=Constant((0,0))

F1 = inner(sigma(u), epsilon(v))*x[0]*dx  +v[0]*lmbda*nabla_div(u)*dx + v[0]*(lmbda+2*mu)*u[0]/x[0]*dx - alpha*p*nabla_div(v)*x[0]*dx -alpha*p*v[0]*dx\
    -inner(T, v)*x[0]*ds(subdomain_id=1, domain=mesh, subdomain_data=contorno)


F2 = dt*inner(nabla_grad(q),nabla_grad(K*p))*x[0]*dx\
+ alpha*divu_t*q*dx +Constant((s_coef))*(dp_t)*q*x[0]*dx\
-dt*(inner(Constant((0,0)),n))*q*ds(subdomain_id=(2,4),domain=mesh, subdomain_data=contorno) - dt*flo*n[0]*q*ds(subdomain_id=1,domain=mesh, subdomain_data=contorno) 

L_momentum =lhs(F1) 
R_momentum =rhs(F1)
L_mass=lhs(F2)
R_mass=rhs(F2)
L=L_momentum+L_mass
R_S=R_momentum+R_mass
snaps=80


Betta= 1#alpha/(alpha**2 + s_coef*(1/B_m+mu/3))
q_ = Caudal/(2*math.pi*H*K)
M=10
r = np.linspace(0.01, R_, snaps)
def pbar(r,s,f_const):
    ro=r/R_
    sigmma = (s * R_ ** 2) / cv_dot
    I_cero_ro = besseli(0, ro * (sigmma) ** (1 / 2))
    K_cero_ro = besselk(0, ro * (sigmma) ** (1 / 2))
    I_cero = besseli(0, (sigmma) ** (1 / 2))
    K_cero = besselk(0, (sigmma) ** (1 / 2))
    term_1 = (Betta * mu * f_const / q_) * (1 - I_cero_ro / I_cero)
    term_2 = (1 / sigmma) *( K_cero_ro-(K_cero / I_cero) * I_cero_ro)
    v_ray = (term_1 - term_2 )
    return v_ray

def f_coef(r,s,f_const):
    ro=r/R_
    sigmma = (s * R_ ** 2)/cv_dot
    part_1 = (2*alpha/sigmma**2)*(1-1/(besseli(0,np.sqrt(sigmma))))
    part_2 = 3/(B_m*mu) + alpha*Betta*(1-(2*besseli(1,np.sqrt(sigmma)))/(np.sqrt(sigmma)*besseli(0,np.sqrt(sigmma))))
    return (part_1/part_2)*q_/mu
def rad_displacement(r,s,f_const):
    ro = r/R_
    sigmma=(s * R_ ** 2) / cv_dot
    I_uno_ro = besseli(1, ro * (sigmma) ** (1 / 2))
    K_uno_ro = besselk(1, ro * (sigmma) ** (1 / 2))
    I_cero = besseli(0, (sigmma) ** (1 / 2))
    K_cero = besselk(0, (sigmma) ** (1 / 2))
    part_1_1 = ro/4- I_uno_ro/(2*sigmma**(1/2)*I_cero)
    part_1_2 = ro/4*(K_m/(mu)+4/3)
    part_1 = (alpha*Betta*part_1_1-part_1_2)*(mu*f_const/q_)
    part_2 = (alpha/(2*sigmma**(3/2)))*(K_uno_ro - 1/(ro*sigmma ** (1/2)) + (K_cero/I_cero)*I_uno_ro)
    n_ray = part_1 + part_2 

    return n_ray

def talbot(M, t,Fun,r=0,f_const=0):
    v=0
    for k in range(M ):
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
    f_const = talbot(M, t, f_coef)
    for l in range(len(r)):   
        if l==0:
            v_final=np.array([[talbot(M, t, pbar,r=r[0],f_const=f_const)*R_**2/cv_dot,r[l]/R_]])
        else:
            v_final= np.append(v_final,np.array([[talbot(M, t, pbar,r=r[l],f_const=f_const)*R_**2/cv_dot,r[l]/R_]]),axis=0)
    return v_final


def u_analitico(r, M, t):
    f_const = talbot(M, t, f_coef)
    for l in range(len(r)):
        if l == 0:
            u_final= np.array([[-talbot(M, t, rad_displacement, r=r[0], f_const=f_const)*R_**2/cv_dot, r[l]/R_]])
        else:
            u_final = np.append(u_final,np.array([[-talbot(M, t, rad_displacement, r=r[l], f_const=f_const)*R_**2/cv_dot, r[l]/R_]]),axis=0)
    return u_final
def u_analiticoP(r, M, t):
    f_const = talbot(M, t, f_coef)
    u_final=-talbot(M, t, rad_displacement, r=r, f_const=f_const)
    return -u_final*R_**2/cv_dot

t=delta
X_w = Function(W)
bcs=[bc1,bc2,bp1]
f = plt.figure()
f.set_figwidth(10)
f.set_figheight(10)
L2=[]
uplot=[]
urplot=[]
ig, pres = plt.subplots()
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
    

    # if pot % 10== 0:
        # s = sigma(u_)
        
        # cauchy=project(s,TS)

        # o1 = sigma_1(cauchy)
        
        # o2 = sigma_3(cauchy)
        
        # tm=(o1-o2)/2
        # fail = envFalla(o1, o2,theta,C)
        # fs = tm
        # fs=project(fs,Z)
        # flow=-K*grad(p_)
        # flow=project(flow,Z_v)
        # fs.rename(" mean stress", "mean stress") ;vtkfile_fs.write(fs, t)
        # u_.rename("displacement", "displacement") ;vtkfile_u.write(u_, t)
        # flow.rename("flow", "flow") ;vtkfile_flow.write(flow, t)
        # p_.rename("pressure", "pressure"); vtkfile_p.write(p_, t)
    uplot.append([(lmbda+mu)*u_(0.015*R_,0)[0]/(R_*q_),(lmbda+mu)*u_(0.13*R_,0)[0]/(R_*q_),(lmbda+mu)*u_(R_,0)[0]/(R_*q_),tdot])
    #urplot.append([u_analiticoP(0.01*R_, M, t),u_analiticoP(0.1*R_, M, t),u_analiticoP(R_, M, t),tdot])
    z_=0
    if near(t,0.01/(cv_dot/R_**2),dt/2) or near(t,0.1/(cv_dot/R_**2),dt/2) or near(t,1/(cv_dot/R_**2),dt/2) :
        for y in range(snaps):
            pdot=(p_(z_*R_,0)/(q_))
            udot=(-(lmbda+mu)*u_(z_*R_,H)[0]/(R_*q_))
            if y ==0:
                results=np.array([[pdot,udot,z_]])
            else:
                results =np.append(results,np.array([[pdot,udot,z_]]),axis=0)
            z_+=1/snaps
        p_a = p_analitico(r,M,t)
        #u_a = u_analitico(r,M,t)
        #L2.append([np.sum((results[:,0]**2-p_a[:,0])**2),tdot])
        #pressure grap
        line1, =pres.plot(results[:, 2],results[:, 0], "-",color='red',label='FEM')
        pres.set_ylabel("$p/q$",fontsize=15)
        pres.set_xlabel("$R/r$",fontsize=15)
        line2,=pres.plot(p_a[:, 1], p_a[:, 0], "--",color='black',label='Analítica')
        pres.set_xlim(0,1)
        pres.set_ylim(-5,0)
        
        # lines = plt.gca().get_lines()
        # l1=lines[-1]
        # labelLine(l1,-2,label=r'$t*=${}'.format(round(tdot,2)),ha='left',va='bottom',align = True)
        # if near(t,0.01/(cv_dot/R_**2),dt/2) :
        #     pres.legend(handler_map={line1: HandlerLine2D(numpoints=4)},loc= 'upper center')
        #displament grap
        # u_a = u_analitico(r,M,t)
        # #pressure grap
        # line1, =disp.plot(results[:, 2],results[:, 1], "-",color='red',label='FEM')
        # disp.set_ylabel("$p/q$",fontsize=20)
        # line2,=disp.plot(u_a[:, 1], u_a[:, 0], "--",color='black',label='Analítica')
        # disp.set_xlim(0,1)
        # disp.set_ylim(0,0.25)
    print('step', pot, 'of', steps,'time*:',tdot)

    t += delta
plt.grid(True,color='k',which="both",alpha=0.3, linestyle='-', linewidth=0.5)
plt.savefig('RESULTADOS_ACUIFERO/presion_aquifero_dt%s_grosse_%s.png'%(round(dtdot,5),scheme),dpi=300)
plt.close()
f = plt.figure()
f.set_figwidth(15.5)
f.set_figheight(7.5)
ig, ax = plt.subplots()

ax.set_ylim(0,0.25)
ax.set_xlim(0.01,10)
uplot=np.array(uplot)
urplot=np.array(urplot)
np.savetxt('RESULTADOS_ACUIFERO/disp%s_grosse_%s.out'%(round(dtdot,5),scheme), (uplot)) 
line1, = ax.plot(uplot[:, 3], -uplot[:,0], "--", color='black')
line2, = ax.plot(uplot[:,3],  -uplot[:,1], "--", color='black')
line3, = ax.plot(uplot[:,3],  -uplot[:,2], "--", color='black',label='FEM')
line4, = ax.plot([-1,-2], [-5,-6], "-", color='red',label='Analitical')
ax.legend(handler_map={line1: HandlerLine2D(numpoints=4)},loc= 'upper left')
ax.set_xscale('log')
plt.xlabel("$t*$ ",fontsize=20)
plt.ylabel("$u*_{z}$",fontsize=20)
plt.grid(True,color='k',which="both",alpha=0.3, linestyle='-', linewidth=0.5)
plt.savefig('RESULTADOS_ACUIFERO/disp_dt%s_grosse_%s.png'%(round(dtdot,5),scheme),dpi=300,)
plt.close()
# u_final=u_analitico(r,M,tdot)
# print(u_final)

# data_001 = np.concatenate(data_001)
# data_05 = np.concatenate(data_05)
# data_1 = np.concatenate(data_1)

# # Trazar la gráfica con las tres líneas
# plt.plot(data_001[:, 1], data_001[:, 0], label='t*=0.01')
# plt.plot(data_05[:, 1], data_05[:, 0], label='t*=0.2')
# plt.plot(data_1[:, 1], data_1[:, 0], label='t*=1')
# plt.xlabel("$r/R$")
# plt.ylabel("$u_z$")
# plt.legend()
# plt.grid(True, color='k', which="both", alpha=0.3, linestyle='-', linewidth=0.5)
# plt.savefig('RESULTADOS_ACUIFERO/desp_aquifero_dt%s_grosse_%s.png'%(round(dtdot,5),scheme),dpi=300)
# plt.close()
# print(radio_001)

# x_values001 = [item[0] for item in radio_001]
# y_values001 = [item[1] for item in radio_001]
# x_values05 = [item[0] for item in radio_05]
# y_values05 = [item[1] for item in radio_05]
# x_values1 = [item[0] for item in radio_1]
# y_values1 = [item[1] for item in radio_1]
# # Trazar la gráfica con las tres líneas
# plt.plot(x_values001, y_values001, label='t*=0.01')
# plt.plot(x_values05, y_values05, label='t*=0.2')
# plt.plot(x_values1, y_values1, label='t*=1')
# plt.xlabel("$t*$")
# plt.ylabel("$u_z$")
# plt.legend()
# plt.grid(True, color='k', which="both", alpha=0.3, linestyle='-', linewidth=0.5)
# plt.savefig('RESULTADOS_ACUIFERO/desp(pvst)_aquifero_dt%s_grosse_%s.png'%(round(dtdot,5),scheme),dpi=300)
# plt.close()