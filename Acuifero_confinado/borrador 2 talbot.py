import math
from mpmath import mpf, besseli, besselk
import numpy as np
import matplotlib.pyplot as plt

steps=1000
R_=300
H=10

# B_s=1E-11
# B_m=1E-10
B_f=4.4E-10
Poro=0.3
nu=1/8#coeficiente de poisson  

E=30000000# #modulo elasticidad
B_m=(E/(3*(1-2*nu)))**(-1)
B_s=B_m/10
alpha=(1-B_s/B_m) #biotcoef
mu = E/(2*(1+nu))#coeficientes de Lame
print("relacion K/G= ",1/(B_m*mu))
lmbda = E*nu/((1+nu)*(1-2*nu))#coeficientes de Lame
mv=1/((1/B_m)+((4/3)*mu))
rango = 100
K=4.5E-4
s_coef=(alpha-Poro)*B_s +Poro*B_f
cv_dot =(K*(1/B_m+mu/3))/(alpha**2+s_coef*(1/B_m+mu/3))
Betta= alpha/(alpha**2 + s_coef*(1/B_m+mu/3))
Caudal=200
q = Caudal/(2*math.pi*R_*H*K)
def pbar_simple (r,s):
    ro = r/R_
    sigmma = (s * R_ ** 2) / cv_dot
    I_cero = besseli(0, (sigmma) ** (1 / 2))
    K_cero = besselk(0, (sigmma) ** (1 / 2))
    I_cero_ro = besseli(0, ro * (sigmma) ** (1 / 2))
    K_cero_ro = besselk(0, ro * (sigmma) ** (1 / 2))
    return -(1 / sigmma) *( K_cero_ro-(K_cero / I_cero) * I_cero_ro)

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

def f_coef(ro,s,f_const):
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
    v_final = []
    f_const = talbot(M, t, f_coef)
    print(f_const)
    print("p analico init")
    for l in range(len(r)):   
        v = 0
        v_final.append([cv_dot*talbot(M, t, pbar,r=r[l],f_const=f_const)/R_**2])
    return np.array(v_final)
def p_analitico_simple(r,t):
    v_final = []

    for l in range(len(r)):   
        v = 0
        v_final.append([t/2*pbar_simple(r[l],t/2)])
    return np.array(v_final)
t=10*R_**2/cv_dot
M=10
r = np.linspace(0.01, R_, 50)

def u_analitico(r,M,t):
    print("u analico init")
    v_final = []
    f_const = talbot(M, t, f_coef)
    for l in range(len(r)):   
        v = 0
        v_final.append([talbot(M, t, rad_displacement,r[l],f_const)])
    return np.array(v_final)
    print("u analico end")
    return v_final
plt.figure(figsize=(10.33,5.25))
for i in range(1):
    print("tiempo adimensional: ",cv_dot*t/R_**2)
    resultado = R_**2*p_analitico(r, M, t)/cv_dot
    plt.plot(r/R_,resultado)
    #plt.plot(r/R_,p_analitico_simple(r, t))
print(R_**2/cv_dot)
#plt.ylim((0,0.25))
#plt.xlim((0,1))
plt.plot(r/R_,np.log(r/R_))


